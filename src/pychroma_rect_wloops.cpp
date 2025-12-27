#include "pychroma_api.h"
#include "pychroma_state.h"
#include "pychroma_util.h"

#include "util/ft/time_slice_set.h"

#include <hdf5.h>
#include <vector>

extern "C" {

int pychroma_rect_wloops(int t_dir, int z_dir, int L_max, int R_max, const char* out_path) {
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
  if (out_path == nullptr) {
    set_error("Output path is null.");
    return 1;
  }
  if (t_dir < 0 || t_dir >= Nd || z_dir < 0 || z_dir >= Nd || t_dir == z_dir) {
    set_error("Invalid t_dir/z_dir.");
    return 1;
  }
  if (L_max <= 0 || R_max <= 0) {
    set_error("L_max and R_max must be positive.");
    return 1;
  }

  const int Nt = Layout::lattSize()[t_dir];
  if (L_max >= Nt) {
    set_error("L_max must be less than the time extent.");
    return 1;
  }

  int transverse[2];
  int idx = 0;
  for (int d = 0; d < Nd; ++d) {
    if (d != t_dir && d != z_dir) {
      transverse[idx++] = d;
    }
  }
  if (idx != 2) {
    set_error("Failed to determine transverse directions.");
    return 1;
  }
  const int x_dir = transverse[0];
  const int y_dir = transverse[1];

  const double spatial_vol = double(Layout::vol()) / double(Nt);
  const double norm = 1.0 / (2.0 * spatial_vol * double(Nc));

  TimeSliceSet tset(t_dir);
  const Set& tset_set = tset.getSet();

  std::vector<LatticeColorMatrix> z_paths = path_products_dir(g_u, z_dir, L_max);
  std::vector<LatticeColorMatrix> x_paths = path_products_dir(g_u, x_dir, R_max);
  std::vector<LatticeColorMatrix> y_paths = path_products_dir(g_u, y_dir, R_max);

  struct H5Complex {
    double re;
    double im;
  };
  const size_t total = size_t(Nt) * size_t(L_max) * size_t(R_max);
  std::vector<H5Complex> data(total);

  for (int L = 1; L <= L_max; ++L) {
    const LatticeColorMatrix& U_z = z_paths[L];
    for (int R = 1; R <= R_max; ++R) {
      multi1d<DComplex> sum_t(Nt);
      for (int t = 0; t < Nt; ++t) {
        sum_t[t] = cmplx(Real(0), Real(0));
      }

      const LatticeColorMatrix& U_x = x_paths[R];
      const LatticeColorMatrix& U_y = y_paths[R];
      const LatticeColorMatrix U_x_shift = shift_n(U_x, z_dir, L);
      const LatticeColorMatrix U_y_shift = shift_n(U_y, z_dir, L);
      const LatticeColorMatrix U_z_shift_x = shift_n(U_z, x_dir, R);
      const LatticeColorMatrix U_z_shift_y = shift_n(U_z, y_dir, R);

      {
        LatticeColorMatrix loop = U_z * U_x_shift * adj(U_z_shift_x) * adj(U_x);
        LatticeComplex tr = trace(loop);
        multi1d<DComplex> sum_dir = sumMulti(tr, tset_set);
        for (int t = 0; t < Nt; ++t) {
          sum_t[t] += sum_dir[t];
        }
      }
      {
        LatticeColorMatrix loop = U_z * U_y_shift * adj(U_z_shift_y) * adj(U_y);
        LatticeComplex tr = trace(loop);
        multi1d<DComplex> sum_dir = sumMulti(tr, tset_set);
        for (int t = 0; t < Nt; ++t) {
          sum_t[t] += sum_dir[t];
        }
      }

      for (int t = 0; t < Nt; ++t) {
        DComplex val = sum_t[t] * Real(norm);
        const size_t idx_out = size_t(t) + size_t(Nt) * (size_t(L - 1) * size_t(R_max) + size_t(R - 1));
        data[idx_out].re = toDouble(real(val));
        data[idx_out].im = toDouble(imag(val));
      }
    }
  }

  hid_t file = H5Fcreate(out_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file < 0) {
    set_error("Failed to create HDF5 file.");
    return 1;
  }

  hid_t complex_type = H5Tcreate(H5T_COMPOUND, sizeof(H5Complex));
  H5Tinsert(complex_type, "r", HOFFSET(H5Complex, re), H5T_NATIVE_DOUBLE);
  H5Tinsert(complex_type, "i", HOFFSET(H5Complex, im), H5T_NATIVE_DOUBLE);

  hsize_t dims[3] = {static_cast<hsize_t>(Nt), static_cast<hsize_t>(L_max), static_cast<hsize_t>(R_max)};
  hid_t space = H5Screate_simple(3, dims, nullptr);
  hid_t dset = H5Dcreate2(file, "wilson_loops", complex_type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dset < 0) {
    H5Sclose(space);
    H5Tclose(complex_type);
    H5Fclose(file);
    set_error("Failed to create HDF5 dataset.");
    return 1;
  }

  if (H5Dwrite(dset, complex_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()) < 0) {
    H5Dclose(dset);
    H5Sclose(space);
    H5Tclose(complex_type);
    H5Fclose(file);
    set_error("Failed to write HDF5 dataset.");
    return 1;
  }

  auto write_attr_int = [&](const char* name, int value) {
    hid_t attr_space = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate2(file, name, H5T_NATIVE_INT, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr >= 0) {
      H5Awrite(attr, H5T_NATIVE_INT, &value);
      H5Aclose(attr);
    }
    H5Sclose(attr_space);
  };

  write_attr_int("t_dir", t_dir);
  write_attr_int("z_dir", z_dir);
  write_attr_int("L_max", L_max);
  write_attr_int("R_max", R_max);
  write_attr_int("Nc", Nc);

  H5Dclose(dset);
  H5Sclose(space);
  H5Tclose(complex_type);
  H5Fclose(file);
  return 0;
}

} // extern "C"
