#include "pychroma_api.h"
#include "pychroma_state.h"

#include "util/gauge/gauge_init_aggregate.h"
#include "util/gauge/rgauge.h"
#include "io/xml_group_reader.h"
#include "meas/inline/inline.h"
#include "meas/inline/io/default_gauge_field.h"

#include <cstdio>
#include <cstdlib>
#include <sstream>

extern "C" {

int pychroma_initialize() {
  if (g_initialized) {
    return 0;
  }
  const char* pool_env = std::getenv("PYCHROMA_POOL_GB");
  const char* pool_gb = pool_env ? pool_env : "0.25";
  const char* out_env = std::getenv("PYCHROMA_OUT");
  int argc = out_env ? 5 : 3;
  static char arg0[] = "py_chroma";
  static char pool_arg[] = "-poolsize";
  static char pool_val[32];
  std::snprintf(pool_val, sizeof(pool_val), "%s", pool_gb);
  static char out_arg[] = "-o";
  static char out_val[1024];
  if (out_env) {
    std::snprintf(out_val, sizeof(out_val), "%s", out_env);
  }
  static char* argv_list[] = {arg0, pool_arg, pool_val, out_arg, out_val, nullptr};
  char** argv = argv_list;
  Chroma::initialize(&argc, &argv);
  g_initialized = true;
  return 0;
}

int pychroma_finalize() {
  if (!g_initialized) {
    return 0;
  }
  Chroma::finalize();
  g_initialized = false;
  g_layout_ready = false;
  g_gauge_ready = false;
  g_output_started = false;
  g_u.resize(0);
  g_input_xml.clear();
  g_gauge_xml.clear();
  return 0;
}

int pychroma_register_inline() {
  if (guard_initialized() != 0) {
    return 1;
  }
  InlineAggregateEnv::registerAll();
  GaugeInitEnv::registerAll();
  return 0;
}

int pychroma_set_lattice(int nd, const int* nrow) {
  if (guard_initialized() != 0) {
    return 1;
  }
  if (nd <= 0 || nrow == nullptr) {
    set_error("Invalid lattice dimensions.");
    return 1;
  }
  multi1d<int> sizes(nd);
  for (int i = 0; i < nd; ++i) {
    sizes[i] = nrow[i];
  }
  Layout::setLattSize(sizes);
  Layout::create();
  g_layout_ready = true;
  return 0;
}

int pychroma_set_rng_seed(unsigned long seed) {
  if (guard_initialized() != 0) {
    return 1;
  }
  QDP::RNG::setrn(seed);
  return 0;
}

int pychroma_set_rng_seed_xml(const char* rng_xml, const char* rng_path) {
  if (guard_initialized() != 0) {
    return 1;
  }
  if (rng_xml == nullptr || rng_path == nullptr) {
    set_error("Invalid RNG XML inputs.");
    return 1;
  }

  try {
    std::istringstream xml_stream(rng_xml);
    XMLReader rng_reader(xml_stream);
    QDP::Seed seed;
    read(rng_reader, rng_path, seed);
    QDP::RNG::setrn(seed);
  } catch (const std::string& e) {
    set_error(std::string("RNG init failed: ") + e);
    return 1;
  } catch (const std::exception& e) {
    set_error(std::string("RNG init failed: ") + e.what());
    return 1;
  } catch (...) {
    set_error("RNG init failed: unknown error.");
    return 1;
  }
  return 0;
}

int pychroma_set_gauge(const char* cfg_id, const char* cfg_xml, const char* cfg_path) {
  if (guard_initialized() != 0) {
    return 1;
  }
  if (!g_layout_ready) {
    set_error("Layout not initialized. Call pychroma_set_lattice() first.");
    return 1;
  }
  if (cfg_id == nullptr || cfg_xml == nullptr || cfg_path == nullptr) {
    set_error("Invalid gauge configuration inputs.");
    return 1;
  }

  try {
    std::istringstream xml_stream(cfg_xml);
    XMLReader cfgtop(xml_stream);

    Handle<GaugeInit> gaugeInit(TheGaugeInitFactory::Instance().createObject(
        std::string(cfg_id), cfgtop, std::string(cfg_path)));

    XMLReader gauge_file_xml;
    XMLReader gauge_xml;
    g_u.resize(Nd);
    (*gaugeInit)(gauge_file_xml, gauge_xml, g_u);

    XMLBufferWriter config_xml;
    config_xml << gauge_xml;
    g_gauge_xml = config_xml.str();

    InlineDefaultGaugeField::reset();
    InlineDefaultGaugeField::set(g_u, config_xml);
    g_gauge_ready = true;
  } catch (const std::string& e) {
    set_error(std::string("Gauge init failed: ") + e);
    return 1;
  } catch (const std::exception& e) {
    set_error(std::string("Gauge init failed: ") + e.what());
    return 1;
  } catch (...) {
    set_error("Gauge init failed: unknown error.");
    return 1;
  }

  return 0;
}

int pychroma_set_input_xml(const char* input_xml) {
  if (input_xml == nullptr) {
    set_error("Input XML is null.");
    return 1;
  }
  g_input_xml = input_xml;
  return 0;
}

int pychroma_random_gauge_transform() {
  if (guard_initialized() != 0) {
    return 1;
  }
  if (!g_gauge_ready) {
    set_error("Gauge field not set. Call pychroma_set_gauge() first.");
    return 1;
  }

  try {
    rgauge(g_u);

    XMLBufferWriter config_xml;
    if (!g_gauge_xml.empty()) {
      std::istringstream gauge_stream(g_gauge_xml);
      XMLReader gauge_reader(gauge_stream);
      config_xml << gauge_reader;
    }
    InlineDefaultGaugeField::reset();
    InlineDefaultGaugeField::set(g_u, config_xml);
  } catch (const std::string& e) {
    set_error(std::string("Random gauge transform failed: ") + e);
    return 1;
  } catch (const std::exception& e) {
    set_error(std::string("Random gauge transform failed: ") + e.what());
    return 1;
  } catch (...) {
    set_error("Random gauge transform failed: unknown error.");
    return 1;
  }

  return 0;
}

const char* pychroma_last_error() {
  return g_last_error.c_str();
}

} // extern "C"
