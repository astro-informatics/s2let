// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <complex.h>
#include <ssht.h>
#include <stdlib.h>
#include <math.h>

/*!
 * Interface to SSHT (required by the Java interface to S2LET)
 * Output : MW complex signal
 *
 * \param[out]  f Output MW map.
 * \param[in]  flm Spherical harmonic coefficients.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_mw_alm2map(complex double* f, const complex double* flm, int L, int spin) {
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;
  ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);
}

/*!
 * Interface to SSHT (required by the Java interface to S2LET)
 * Output : spherical harmonics
 *
 * \param[out]  flm Spherical harmonic coefficients.
 * \param[in]  f Output MW map.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_mw_map2alm(complex double* flm, const complex double* f, int L, int spin) {
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;
  ssht_core_mw_forward_sov_conv_sym(flm, f, L, spin, dl_method, verbosity);
}

/*!
 * Interface to SSHT (required by the Java interface to S2LET)
 * Output : MW real signal
 *
 * \param[out]  f Output MW map.
 * \param[in]  flm Spherical harmonic coefficients.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_mw_alm2map_real(double* f, const complex double* flm, int L) {
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;
  ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);
}

/*!
 * Interface to SSHT (required by the Java interface to S2LET)
 * Output : spherical harmonics (corresponding to real signal)
 *
 * \param[out]  flm Spherical harmonic coefficients.
 * \param[in]  f Output MW map.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_mw_map2alm_real(complex double* flm, const double* f, int L) {
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;
  ssht_core_mw_forward_sov_conv_sym_real(flm, f, L, dl_method, verbosity);
}

/*!
 * Computes power of a complex MW signal
 */
double s2let_mw_power(complex double *f, int L){
  complex double *flm;
  s2let_allocate_lm(&flm, L);
  s2let_mw_map2alm(flm, f, L, 0);
  double res = s2let_lm_power(flm, L);
  free(flm);
  return res;
}

/*!
 * Computes power of a real MW signal
 */
double s2let_mw_power_real(double *f, int L){
  complex double *flm;
  s2let_allocate_lm(&flm, L);
  s2let_mw_map2alm_real(flm, f, L);
  double res = s2let_lm_power(flm, L);
  free(flm);
  return res;
}
