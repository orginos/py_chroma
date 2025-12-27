#include "pychroma_api.h"
#include "pychroma_state.h"

#include "meas/inline/inline.h"
#include "meas/inline/io/default_gauge_field.h"

#include <sstream>

extern "C" {

int pychroma_run_inline_xml(const char* inline_xml) {
  if (guard_initialized() != 0) {
    return 1;
  }
  if (!g_layout_ready) {
    set_error("Layout not initialized. Call pychroma_set_lattice() first.");
    return 1;
  }
  if (!g_gauge_ready) {
    set_error("Gauge field not set. Call pychroma_set_gauge() first.");
    return 1;
  }
  if (inline_xml == nullptr) {
    set_error("Inline XML string is null.");
    return 1;
  }

  try {
    std::istringstream meas_stream(inline_xml);
    XMLReader meas_xml(meas_stream);
    multi1d<Handle<AbsInlineMeasurement> > the_measurements;
    read(meas_xml, "/InlineMeasurements", the_measurements);

    XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
    if (!g_output_started) {
      push(xml_out, "chroma");
      if (!g_input_xml.empty()) {
        std::istringstream input_stream(g_input_xml);
        XMLReader input_reader(input_stream);
        write(xml_out, "Input", input_reader);
      }
      proginfo(xml_out);
      QDP::Seed seed;
      QDP::RNG::savern(seed);
      write(xml_out, "RNG", seed);
      if (!g_gauge_xml.empty()) {
        std::istringstream gauge_stream(g_gauge_xml);
        XMLReader gauge_reader(gauge_stream);
        write(xml_out, "Config_info", gauge_reader);
      }
      MesPlq(xml_out, "Observables", g_u);
      g_output_started = true;
    }

    push(xml_out, "InlineObservables");
    xml_out.flush();

    unsigned long cur_update = 0;
    for (int m = 0; m < the_measurements.size(); ++m) {
      AbsInlineMeasurement& meas = *(the_measurements[m]);
      if (cur_update % meas.getFrequency() == 0) {
        push(xml_out, "elem");
        meas(cur_update, xml_out);
        pop(xml_out);
        xml_out.flush();
      }
    }

    pop(xml_out);
    if (g_output_started) {
      pop(xml_out);
    }
  } catch (const std::string& e) {
    set_error(std::string("Inline run failed: ") + e);
    return 1;
  } catch (const std::exception& e) {
    set_error(std::string("Inline run failed: ") + e.what());
    return 1;
  } catch (...) {
    set_error("Inline run failed: unknown error.");
    return 1;
  }

  return 0;
}

int pychroma_run_plaquette(unsigned long update_no, unsigned long frequency) {
  if (guard_initialized() != 0) {
    return 1;
  }
  if (!g_layout_ready) {
    set_error("Layout not initialized. Call pychroma_set_lattice() first.");
    return 1;
  }
  if (!g_gauge_ready) {
    set_error("Gauge field not set. Call pychroma_set_gauge() first.");
    return 1;
  }
  if (frequency == 0) {
    frequency = 1;
  }

  try {
    InlinePlaquetteEnv::Params params;
    params.frequency = frequency;
    params.named_obj.gauge_id = InlineDefaultGaugeField::getId();
    params.xml_file = "";

    InlinePlaquetteEnv::InlineMeas meas(params);

    XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
    push(xml_out, "InlineObservables");
    push(xml_out, "elem");
    meas(update_no, xml_out);
    pop(xml_out);
    pop(xml_out);
    xml_out.flush();
  } catch (const std::string& e) {
    set_error(std::string("Plaquette failed: ") + e);
    return 1;
  } catch (const std::exception& e) {
    set_error(std::string("Plaquette failed: ") + e.what());
    return 1;
  } catch (...) {
    set_error("Plaquette failed: unknown error.");
    return 1;
  }

  return 0;
}

} // extern "C"
