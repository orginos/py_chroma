#include "pychroma_util.h"

LatticeColorMatrix shift_n(const LatticeColorMatrix& in, int dir, int n) {
  LatticeColorMatrix out = in;
  for (int i = 0; i < n; ++i) {
    out = shift(out, FORWARD, dir);
  }
  return out;
}

std::vector<LatticeColorMatrix> path_products_dir(const multi1d<LatticeColorMatrix>& u, int dir, int max_len) {
  std::vector<LatticeColorMatrix> paths(max_len + 1);
  paths[0] = 1;
  if (max_len == 0) {
    return paths;
  }
  LatticeColorMatrix link = u[dir];
  paths[1] = link;
  for (int len = 2; len <= max_len; ++len) {
    link = shift(link, FORWARD, dir);
    paths[len] = paths[len - 1] * link;
  }
  return paths;
}
