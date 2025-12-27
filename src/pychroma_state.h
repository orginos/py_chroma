#ifndef PYCHROMA_STATE_H
#define PYCHROMA_STATE_H

#include "chroma.h"

#include <string>

extern std::string g_last_error;
extern bool g_initialized;
extern bool g_layout_ready;
extern bool g_gauge_ready;
extern bool g_output_started;
extern multi1d<LatticeColorMatrix> g_u;
extern std::string g_input_xml;
extern std::string g_gauge_xml;

void set_error(const std::string& msg);
int guard_initialized();

#endif
