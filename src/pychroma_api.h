#ifndef PYCHROMA_API_H
#define PYCHROMA_API_H

#ifdef __cplusplus
extern "C" {
#endif

int pychroma_initialize();
int pychroma_finalize();
int pychroma_register_inline();
int pychroma_set_lattice(int nd, const int* nrow);
int pychroma_set_rng_seed(unsigned long seed);
int pychroma_set_rng_seed_xml(const char* rng_xml, const char* rng_path);
int pychroma_set_gauge(const char* cfg_id, const char* cfg_xml, const char* cfg_path);
int pychroma_random_gauge_transform();
int pychroma_set_input_xml(const char* input_xml);
int pychroma_run_inline_xml(const char* inline_xml);
int pychroma_run_plaquette(unsigned long update_no, unsigned long frequency);
int pychroma_rect_wloops(int t_dir, int z_dir, int L_max, int R_max, const char* out_path);
int pychroma_run_hmc_xml(const char* params_xml);
int pychroma_run_smd_xml(const char* params_xml);
const char* pychroma_last_error();

#ifdef __cplusplus
}
#endif

#endif
