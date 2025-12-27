#ifndef PYCHROMA_UTIL_H
#define PYCHROMA_UTIL_H

#include "chroma.h"

#include <vector>

LatticeColorMatrix shift_n(const LatticeColorMatrix& in, int dir, int n);
std::vector<LatticeColorMatrix> path_products_dir(const multi1d<LatticeColorMatrix>& u, int dir, int max_len);

#endif
