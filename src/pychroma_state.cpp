#include "pychroma_state.h"

std::string g_last_error;
bool g_initialized = false;
bool g_layout_ready = false;
bool g_gauge_ready = false;
bool g_output_started = false;
multi1d<LatticeColorMatrix> g_u;
std::string g_input_xml;
std::string g_gauge_xml;

void set_error(const std::string& msg) {
  g_last_error = msg;
}

int guard_initialized() {
  if (!g_initialized) {
    set_error("Chroma not initialized. Call pychroma_initialize() first.");
    return 1;
  }
  return 0;
}
