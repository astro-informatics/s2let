#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>

#include "s2let.h"

#include <cmocka.h>

const int seed = 1;

_Bool min(int a, int b) { return a > b ? b : a; }

void test_transform_axisym_vs_directional_mw(void **state) {
  s2let_parameters_t parameters = *(s2let_parameters_t *)*state;
  parameters.upsample = 1;

  int J = s2let_j_max(&parameters);

  double wav_error, scal_error;

  complex double *f, *flm;
  s2let_allocate_lm(&flm, parameters.L);
  s2let_allocate_mw(&f, parameters.L);
  s2let_lm_random_flm(flm, parameters.L, parameters.spin, seed);

  ssht_core_mw_inverse_sov_sym(
      f,
      flm,
      parameters.L,
      parameters.spin,
      parameters.dl_method,
      parameters.verbosity);

  complex double *f_wav_axisym, *f_scal_axisym, *f_wav_dir, *f_scal_dir;
  s2let_transform_axisym_allocate_mw_f_wav(&f_wav_axisym, &f_scal_axisym, &parameters);
  s2let_allocate_f_wav(&f_wav_dir, &f_scal_dir, &parameters);

  s2let_transform_axisym_wav_analysis_mw(f_wav_axisym, f_scal_axisym, f, &parameters);
  s2let_analysis_px2wav(f_wav_dir, f_scal_dir, f, &parameters);

  // Compute the maximum absolute error in the computed wavelet transform
  for (int i = 0;
       i < (J - parameters.J_min + 1) * parameters.L * (2 * parameters.L - 1);
       i += 1) {
    assert_float_equal(creal(f_wav_axisym[i]), creal(f_wav_dir[i]), 1e-12);
    assert_float_equal(cimag(f_wav_axisym[i]), cimag(f_wav_dir[i]), 1e-12);
  }
  for (int i = 0; i < parameters.L * (2 * parameters.L - 1); i += 1) {
    assert_float_equal(creal(f_scal_axisym[i]), creal(f_scal_dir[i]), 1e-12);
    assert_float_equal(cimag(f_scal_axisym[i]), cimag(f_scal_dir[i]), 1e-12);
  }
}

void test_transform_axisym_vs_directional_mw_multires(void **state) {
  const s2let_parameters_t parameters = *(s2let_parameters_t *)*state;

  const int J = s2let_j_max(&parameters);
  double wav_error, scal_error;
  complex double *f, *flm;
  s2let_allocate_lm(&flm, parameters.L);
  s2let_allocate_mw(&f, parameters.L);
  s2let_lm_random_flm(flm, parameters.L, parameters.spin, seed);

  ssht_core_mw_inverse_sov_sym(
      f,
      flm,
      parameters.L,
      parameters.spin,
      parameters.dl_method,
      parameters.verbosity);

  complex double *f_wav_axisym, *f_scal_axisym, *f_wav_dir, *f_scal_dir;
  s2let_transform_axisym_allocate_mw_f_wav_multires(
      &f_wav_axisym, &f_scal_axisym, &parameters);
  s2let_allocate_f_wav(&f_wav_dir, &f_scal_dir, &parameters);

  s2let_transform_axisym_wav_analysis_mw_multires(
      f_wav_axisym, f_scal_axisym, f, &parameters);
  s2let_analysis_px2wav(f_wav_dir, f_scal_dir, f, &parameters);

  int samples = 0;
  for (int j = parameters.J_min; j <= J; ++j) {
    const int bandlim = min(s2let_bandlimit(j, &parameters), parameters.L);
    samples += bandlim * (2 * bandlim - 1);
  }

  // Compute the maximum absolute error in the computed wavelet transform
  for (int i = 0; i < samples; i += 1) {
    assert_float_equal(creal(f_wav_axisym[i]), creal(f_wav_dir[i]), 1e-12);
    assert_float_equal(cimag(f_wav_axisym[i]), cimag(f_wav_dir[i]), 1e-12);
  }
  const int bandlimit =
      min(s2let_bandlimit(parameters.J_min - 1, &parameters), parameters.L);
  for (int i = 0; i < bandlimit * (2 * bandlimit - 1); i += 1) {
    assert_float_equal(creal(f_scal_axisym[i]), creal(f_scal_dir[i]), 1e-12);
    assert_float_equal(cimag(f_scal_axisym[i]), cimag(f_scal_dir[i]), 1e-12);
  }

  free(flm);
  free(f);
  free(f_wav_axisym);
  free(f_scal_axisym);
  free(f_wav_dir);
  free(f_scal_dir);
}

int main(void) {
  const s2let_parameters_t state = {
      .L = 16,
      .N = 1,
      .B = 2,
      .J_min = 2,
      .spin = 0,
      .original_spin = 0,
      .upsample = 0,
      .dl_method = SSHT_DL_RISBO,
      .reality = 0,
      .sampling_scheme = S2LET_SAMPLING_MW,
  };
  const struct CMUnitTest tests[] = {
      cmocka_unit_test_prestate(
          test_transform_axisym_vs_directional_mw, (void **)&state),
      cmocka_unit_test_prestate(
          test_transform_axisym_vs_directional_mw_multires, (void **)&state),
  };

  return cmocka_run_group_tests(tests, NULL, NULL);
}
