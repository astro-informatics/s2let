#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>

#include "s2let.h"

#include <cmocka.h>

const int seed = 1;

void test_transform_axisym_lm_wav(void **state) {
  const s2let_parameters_t *const parameters = (const s2let_parameters_t *)*state;

  double *wav_lm, *scal_lm;

  // Allocate the wavelet kernels
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, parameters);

  // Compute the wavelet kernels
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, parameters);

  complex double *f_wav_lm, *f_scal_lm, *flm, *flm_rec;
  s2let_allocate_lm(&flm, parameters->L);
  s2let_allocate_lm(&flm_rec, parameters->L);

  // Generate a random spherical harmonic decomposition
  s2let_lm_random_flm(flm, parameters->L, 0, seed);

  // Allocate space for the wavelet scales (their harmonic coefficients)
  s2let_transform_axisym_lm_allocate_f_wav(&f_wav_lm, &f_scal_lm, parameters);

  // Perform the wavelet transform through exact harmonic tiling
  s2let_transform_axisym_lm_wav_analysis(
      f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, parameters);

  // Reconstruct the initial harmonic coefficients from those of the wavelets
  s2let_transform_axisym_lm_wav_synthesis(
      flm_rec, f_wav_lm, f_scal_lm, wav_lm, scal_lm, parameters);

  // Compute the maximum absolute error on the harmonic coefficients
  for (int i = 0; i < parameters->L * parameters->L; i += 1) {
    assert_float_equal(creal(flm[i]), creal(flm_rec[i]), 1e-12);
    assert_float_equal(cimag(flm[i]), cimag(flm_rec[i]), 1e-12);
  }

  free(flm);
  free(flm_rec);
  free(f_wav_lm);
  free(f_scal_lm);
  free(wav_lm);
  free(scal_lm);
}

void test_transform_axisym_lm_wav_multires(void **state) {
  const s2let_parameters_t *const parameters = (const s2let_parameters_t *)*state;

  double *wav_lm, *scal_lm;

  // Allocate the wavelet kernels
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, parameters);

  // Compute the wavelet kernels
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, parameters);

  complex double *f_wav_lm, *f_scal_lm, *flm, *flm_rec;
  s2let_allocate_lm(&flm, parameters->L);
  s2let_allocate_lm(&flm_rec, parameters->L);

  // Generate a random spherical harmonic decomposition
  s2let_lm_random_flm(flm, parameters->L, 0, seed);

  // Allocate space for the wavelet scales (their harmonic coefficients)
  s2let_transform_axisym_lm_allocate_f_wav_multires(&f_wav_lm, &f_scal_lm, parameters);

  // Perform the wavelet transform through exact harmonic tiling
  s2let_transform_axisym_lm_wav_analysis_multires(
      f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, parameters);

  // Reconstruct the initial harmonic coefficients from those of the wavelets
  s2let_transform_axisym_lm_wav_synthesis_multires(
      flm_rec, f_wav_lm, f_scal_lm, wav_lm, scal_lm, parameters);

  // Compute the maximum absolute error on the harmonic coefficients
  for (int i = 0; i < parameters->L * parameters->L; i += 1) {
    assert_float_equal(creal(flm[i]), creal(flm_rec[i]), 1e-12);
    assert_float_equal(cimag(flm[i]), cimag(flm_rec[i]), 1e-12);
  }

  free(flm);
  free(flm_rec);
  free(f_wav_lm);
  free(f_scal_lm);
  free(wav_lm);
  free(scal_lm);
}

void test_wav_transform_wavlm_manual(void **state) {
  s2let_parameters_t *const parameters = (s2let_parameters_t *)*state;

  const int J = s2let_j_max(parameters) - parameters->J_min;

  complex double *psi;
  double *phi;
  s2let_tiling_wavelet_allocate(&psi, &phi, parameters);
  s2let_tiling_wavelet(psi, phi, parameters);

  const int scal_bandlimit =
      parameters->upsample
          ? parameters->L
          : MIN(s2let_bandlimit(parameters->J_min - 1, parameters), parameters->L);

  complex double *wav_l = (complex double *)calloc(
      (J + 1) * parameters->L * parameters->L, sizeof(complex double));
  int *wav_bandlimits = (int *)calloc(J + 1, sizeof(int));
  for (int j = 0; j <= J; ++j) {
    wav_bandlimits[j] =
        parameters->upsample
            ? parameters->L
            : MIN(s2let_bandlimit(parameters->J_min + j, parameters), parameters->L);
    for (int i = 0; i < parameters->L * parameters->L; ++i)
      wav_l[j * parameters->L * parameters->L + i] =
          psi[(parameters->J_min + j) * parameters->L * parameters->L + i];
  }

  complex double *f_wav, *f_scal, *flm, *flm_rec;
  s2let_allocate_lm(&flm, parameters->L);
  s2let_allocate_lm(&flm_rec, parameters->L);

  // Generate a random spherical harmonic decomposition
  s2let_lm_random_flm(flm, parameters->L, parameters->spin, seed);

  // Allocate space for the wavelet scales (their harmonic/Wigner coefficients)
  s2let_allocate_f_wav_manual(
      &f_wav, &f_scal, wav_bandlimits, scal_bandlimit, parameters->N, J, parameters);

  // Perform the wavelet transform through exact harmonic tiling
  s2let_analysis_lm2wav_manual(
      f_wav,
      f_scal,
      flm,
      phi,
      wav_l,
      scal_bandlimit,
      wav_bandlimits,
      J,
      parameters->L,
      parameters->spin,
      parameters->N);

  // Reconstruct the initial harmonic coefficients from those of the wavelets
  s2let_synthesis_wav2lm_manual(
      flm_rec,
      f_wav,
      f_scal,
      phi,
      wav_l,
      scal_bandlimit,
      wav_bandlimits,
      J,
      parameters->L,
      parameters->spin,
      parameters->N);

  for (int i = 0; i < parameters->L * parameters->L; i += 1) {
    assert_float_equal(creal(flm[i]), creal(flm_rec[i]), 1e-12);
    assert_float_equal(cimag(flm[i]), cimag(flm_rec[i]), 1e-12);
  }

  free(flm);
  free(flm_rec);
  free(f_wav);
  free(f_scal);
  free(psi);
  free(phi);
}

void test_wav_transform_harmonic(void **state) {
  s2let_parameters_t *const parameters = (s2let_parameters_t *)*state;
  complex double *psi;
  double *phi;
  s2let_tiling_wavelet_allocate(&psi, &phi, parameters);
  s2let_tiling_wavelet(psi, phi, parameters);

  complex double *f_wav_lmn, *f_scal_lm, *flm, *flm_rec;
  s2let_allocate_lm(&flm, parameters->L);
  s2let_allocate_lm(&flm_rec, parameters->L);
  s2let_lm_random_flm(flm, parameters->L, parameters->spin, seed);

  s2let_allocate_lmn_f_wav(&f_wav_lmn, &f_scal_lm, parameters);
  s2let_analysis_lm2lmn(f_wav_lmn, f_scal_lm, flm, psi, phi, parameters);
  s2let_synthesis_lmn2lm(flm_rec, f_wav_lmn, f_scal_lm, psi, phi, parameters);

  for (int i = 0; i < parameters->L * parameters->L; i += 1) {
    assert_float_equal(creal(flm[i]), creal(flm_rec[i]), 1e-12);
    assert_float_equal(cimag(flm[i]), cimag(flm_rec[i]), 1e-12);
  }

  free(flm);
  free(flm_rec);
  free(f_wav_lmn);
  free(f_scal_lm);
  free(psi);
  free(phi);
}

void test_transform_axisym_wav(void **state) {
  s2let_parameters_t *const parameters = (s2let_parameters_t *)*state;
  if (parameters->spin != 0)
    skip();

  complex double *f, *f_rec, *flm, *flm_rec;
  s2let_allocate_lm(&flm, parameters->L);
  s2let_allocate_lm(&flm_rec, parameters->L);
  s2let_allocate_mw(&f, parameters->L);
  s2let_allocate_mw(&f_rec, parameters->L);
  s2let_lm_random_flm(flm, parameters->L, 0, seed);

  ssht_core_mw_inverse_sov_sym(
      f,
      flm,
      parameters->L,
      parameters->spin,
      parameters->dl_method,
      parameters->verbosity);

  complex double *f_wav, *f_scal;
  s2let_transform_axisym_allocate_mw_f_wav(&f_wav, &f_scal, parameters);
  s2let_transform_axisym_wav_analysis_mw(f_wav, f_scal, f, parameters);
  s2let_transform_axisym_wav_synthesis_mw(f_rec, f_wav, f_scal, parameters);

  ssht_core_mw_forward_sov_conv_sym(
      flm_rec,
      f_rec,
      parameters->L,
      parameters->spin,
      parameters->dl_method,
      parameters->verbosity);

  for (int i = 0; i < parameters->L * parameters->L; i += 1) {
    assert_float_equal(creal(flm[i]), creal(flm_rec[i]), 1e-12);
    assert_float_equal(cimag(flm[i]), cimag(flm_rec[i]), 1e-12);
  }

  free(f);
  free(f_rec);
  free(f_wav);
  free(f_scal);
}

void test_transform_axisym_wav_real(void **state) {
  s2let_parameters_t parameters = *(s2let_parameters_t *)*state;
  parameters.reality = 1;

  complex double *flm, *flm_rec;
  double *f, *f_rec;
  s2let_allocate_lm(&flm, parameters.L);
  s2let_allocate_lm(&flm_rec, parameters.L);
  s2let_allocate_mw_real(&f, parameters.L);
  s2let_allocate_mw_real(&f_rec, parameters.L);

  s2let_lm_random_flm_real(flm, parameters.L, seed);

  ssht_core_mw_inverse_sov_sym_real(
      f, flm, parameters.L, parameters.dl_method, parameters.verbosity);

  double *f_wav, *f_scal;
  s2let_transform_axisym_allocate_mw_f_wav_real(&f_wav, &f_scal, &parameters);
  s2let_transform_axisym_wav_analysis_mw_real(f_wav, f_scal, f, &parameters);
  s2let_transform_axisym_wav_synthesis_mw_real(f_rec, f_wav, f_scal, &parameters);

  ssht_core_mw_forward_sov_conv_sym_real(
      flm_rec, f_rec, parameters.L, parameters.dl_method, parameters.verbosity);

  for (int i = 0; i < parameters.L * parameters.L; i += 1) {
    assert_float_equal(creal(flm[i]), creal(flm_rec[i]), 1e-12);
    assert_float_equal(cimag(flm[i]), cimag(flm_rec[i]), 1e-12);
  }

  free(f);
  free(f_rec);
  free(f_wav);
  free(f_scal);
}

void test_transform_axisym_wav_multires(void **state) {
  s2let_parameters_t *const parameters = (s2let_parameters_t *)*state;
  if (parameters->spin != 0)
    skip();

  complex double *f, *f_rec, *flm, *flm_rec;
  s2let_allocate_lm(&flm, parameters->L);
  s2let_allocate_lm(&flm_rec, parameters->L);
  s2let_allocate_mw(&f, parameters->L);
  s2let_allocate_mw(&f_rec, parameters->L);
  s2let_lm_random_flm(flm, parameters->L, 0, seed);

  ssht_core_mw_inverse_sov_sym(
      f,
      flm,
      parameters->L,
      parameters->spin,
      parameters->dl_method,
      parameters->verbosity);

  complex double *f_wav, *f_scal;
  s2let_transform_axisym_allocate_mw_f_wav_multires(&f_wav, &f_scal, parameters);
  s2let_transform_axisym_wav_analysis_mw_multires(f_wav, f_scal, f, parameters);
  s2let_transform_axisym_wav_synthesis_mw_multires(f_rec, f_wav, f_scal, parameters);

  ssht_core_mw_forward_sov_conv_sym(
      flm_rec,
      f_rec,
      parameters->L,
      parameters->spin,
      parameters->dl_method,
      parameters->verbosity);

  for (int i = 0; i < parameters->L * parameters->L; i += 1) {
    assert_float_equal(creal(flm[i]), creal(flm_rec[i]), 1e-12);
    assert_float_equal(cimag(flm[i]), cimag(flm_rec[i]), 1e-12);
  }

  free(f);
  free(f_rec);
  free(f_wav);
  free(f_scal);
}

void test_transform_axisym_wav_multires_real(void **state) {
  s2let_parameters_t parameters = *(s2let_parameters_t *)*state;
  parameters.reality = 1;

  complex double *flm, *flm_rec;
  double *f, *f_rec;
  s2let_allocate_lm(&flm, parameters.L);
  s2let_allocate_lm(&flm_rec, parameters.L);
  s2let_allocate_mw_real(&f, parameters.L);
  s2let_allocate_mw_real(&f_rec, parameters.L);

  s2let_lm_random_flm_real(flm, parameters.L, seed);

  ssht_core_mw_inverse_sov_sym_real(
      f, flm, parameters.L, parameters.dl_method, parameters.verbosity);

  double *f_wav, *f_scal;
  s2let_transform_axisym_allocate_mw_f_wav_multires_real(&f_wav, &f_scal, &parameters);
  s2let_transform_axisym_wav_analysis_mw_multires_real(f_wav, f_scal, f, &parameters);
  s2let_transform_axisym_wav_synthesis_mw_multires_real(
      f_rec, f_wav, f_scal, &parameters);

  ssht_core_mw_forward_sov_conv_sym_real(
      flm_rec, f_rec, parameters.L, parameters.dl_method, parameters.verbosity);

  for (int i = 0; i < parameters.L * parameters.L; i += 1) {
    assert_float_equal(creal(flm[i]), creal(flm_rec[i]), 1e-12);
    assert_float_equal(cimag(flm[i]), cimag(flm_rec[i]), 1e-12);
  }

  free(f);
  free(f_rec);
  free(f_wav);
  free(f_scal);
}

void test_wav_transform_mw(void **state) {
  s2let_parameters_t *const parameters = (s2let_parameters_t *)*state;

  complex double *f, *f_rec, *flm, *flm_rec;
  s2let_allocate_lm(&flm, parameters->L);
  s2let_allocate_lm(&flm_rec, parameters->L);
  s2let_allocate_mw(&f, parameters->L);
  s2let_allocate_mw(&f_rec, parameters->L);

  s2let_lm_random_flm(flm, parameters->L, parameters->spin, seed);

  ssht_core_mw_inverse_sov_sym(
      f,
      flm,
      parameters->L,
      parameters->spin,
      parameters->dl_method,
      parameters->verbosity);

  complex double *f_wav, *f_scal;
  s2let_allocate_f_wav(&f_wav, &f_scal, parameters);
  s2let_analysis_px2wav(f_wav, f_scal, f, parameters);
  s2let_synthesis_wav2px(f_rec, f_wav, f_scal, parameters);

  ssht_core_mw_forward_sov_conv_sym(
      flm_rec,
      f_rec,
      parameters->L,
      parameters->spin,
      parameters->dl_method,
      parameters->verbosity);

  for (int i = 0; i < parameters->L * parameters->L; i += 1) {
    assert_float_equal(creal(flm[i]), creal(flm_rec[i]), 1e-12);
    assert_float_equal(cimag(flm[i]), cimag(flm_rec[i]), 1e-12);
  }

  free(f);
  free(f_rec);
  free(flm);
  free(flm_rec);
  free(f_wav);
  free(f_scal);
}

void test_wav_transform_mw_real(void **state) {
  s2let_parameters_t parameters = *(s2let_parameters_t *)*state;
  parameters.reality = 1;
  if (parameters.spin != 0)
    skip();

  double *f, *f_rec;
  complex double *flm, *flm_rec;
  s2let_allocate_lm(&flm, parameters.L);
  s2let_allocate_lm(&flm_rec, parameters.L);
  s2let_allocate_mw_real(&f, parameters.L);
  s2let_allocate_mw_real(&f_rec, parameters.L);

  s2let_lm_random_flm_real(flm, parameters.L, seed);

  ssht_core_mw_inverse_sov_sym_real(
      f, flm, parameters.L, parameters.dl_method, parameters.verbosity);

  double *f_wav, *f_scal;
  s2let_allocate_f_wav_real(&f_wav, &f_scal, &parameters);
  s2let_analysis_px2wav_real(f_wav, f_scal, f, &parameters);
  s2let_synthesis_wav2px_real(f_rec, f_wav, f_scal, &parameters);

  ssht_core_mw_forward_sov_conv_sym_real(
      flm_rec, f_rec, parameters.L, parameters.dl_method, parameters.verbosity);

  for (int i = 0; i < parameters.L * parameters.L; i += 1) {
    assert_float_equal(creal(flm[i]), creal(flm_rec[i]), 1e-12);
    assert_float_equal(cimag(flm[i]), cimag(flm_rec[i]), 1e-12);
  }

  free(f);
  free(f_rec);
  free(flm);
  free(flm_rec);
  free(f_wav);
  free(f_scal);
}

void test_wav_transform_mwss(void **state) {
  skip();
  s2let_parameters_t parameters = *(s2let_parameters_t *)*state;
  parameters.sampling_scheme = S2LET_SAMPLING_MW_SS;

  complex double *f, *f_rec, *flm, *flm_rec;
  s2let_allocate_lm(&flm, parameters.L);
  s2let_allocate_lm(&flm_rec, parameters.L);
  s2let_allocate_mwss(&f, parameters.L);
  s2let_allocate_mwss(&f_rec, parameters.L);

  s2let_lm_random_flm(flm, parameters.L, parameters.spin, seed);

  ssht_core_mw_inverse_sov_sym_ss(
      f,
      flm,
      parameters.L,
      parameters.spin,
      parameters.dl_method,
      parameters.verbosity);

  complex double *f_wav, *f_scal;
  s2let_allocate_f_wav(&f_wav, &f_scal, &parameters);
  s2let_analysis_px2wav(f_wav, f_scal, f, &parameters);
  s2let_synthesis_wav2px(f_rec, f_wav, f_scal, &parameters);

  ssht_core_mw_forward_sov_conv_sym_ss(
      flm_rec,
      f_rec,
      parameters.L,
      parameters.spin,
      parameters.dl_method,
      parameters.verbosity);

  for (int i = 0; i < parameters.L * parameters.L; i += 1) {
    assert_float_equal(creal(flm[i]), creal(flm_rec[i]), 1e-12);
    assert_float_equal(cimag(flm[i]), cimag(flm_rec[i]), 1e-12);
  }

  free(f);
  free(f_rec);
  free(flm);
  free(flm_rec);
  free(f_wav);
  free(f_scal);
}

void test_wav_transform_mwss_real(void **state) {
  s2let_parameters_t parameters = *(s2let_parameters_t *)*state;
  parameters.sampling_scheme = S2LET_SAMPLING_MW_SS;
  if (parameters.spin != 0)
    skip();

  double *f, *f_rec;
  complex double *flm, *flm_rec;
  s2let_allocate_lm(&flm, parameters.L);
  s2let_allocate_lm(&flm_rec, parameters.L);
  s2let_allocate_mwss_real(&f, parameters.L);
  s2let_allocate_mwss_real(&f_rec, parameters.L);

  s2let_lm_random_flm_real(flm, parameters.L, seed);
  ssht_core_mw_inverse_sov_sym_ss_real(
      f, flm, parameters.L, parameters.dl_method, parameters.verbosity);

  double *f_wav, *f_scal;
  s2let_allocate_f_wav_real(&f_wav, &f_scal, &parameters);
  s2let_analysis_px2wav_real(f_wav, f_scal, f, &parameters);
  s2let_synthesis_wav2px_real(f_rec, f_wav, f_scal, &parameters);

  ssht_core_mw_forward_sov_conv_sym_ss_real(
      flm_rec, f_rec, parameters.L, parameters.dl_method, parameters.verbosity);

  for (int i = 0; i < parameters.L * parameters.L; i += 1) {
    assert_float_equal(creal(flm[i]), creal(flm_rec[i]), 1e-12);
    assert_float_equal(cimag(flm[i]), cimag(flm_rec[i]), 1e-12);
  }

  free(f);
  free(f_rec);
  free(flm);
  free(flm_rec);
  free(f_wav);
  free(f_scal);
}

void test_wav_transform_lm2wav(void **state) {
  const s2let_parameters_t parameters = *(s2let_parameters_t *)*state;
  complex double *flm, *flm_rec;
  s2let_allocate_lm(&flm, parameters.L);
  s2let_allocate_lm(&flm_rec, parameters.L);
  s2let_lm_random_flm(flm, parameters.L, parameters.spin, seed);

  complex double *f_wav, *f_scal;
  s2let_allocate_f_wav(&f_wav, &f_scal, &parameters);

  s2let_analysis_lm2wav(f_wav, f_scal, flm, &parameters);
  s2let_synthesis_wav2lm(flm_rec, f_wav, f_scal, &parameters);

  for (int i = 0; i < parameters.L * parameters.L; i += 1) {
    assert_float_equal(creal(flm[i]), creal(flm_rec[i]), 1e-12);
    assert_float_equal(cimag(flm[i]), cimag(flm_rec[i]), 1e-12);
  }

  free(flm);
  free(flm_rec);
  free(f_wav);
  free(f_scal);
}

void test_wav_transform_lm2wav_real(void **state) {
  s2let_parameters_t parameters = *(s2let_parameters_t *)*state;
  parameters.reality = 1;
  if (parameters.spin != 0)
    skip();

  complex double *flm, *flm_rec;
  s2let_allocate_lm(&flm, parameters.L);
  s2let_allocate_lm(&flm_rec, parameters.L);
  s2let_lm_random_flm_real(flm, parameters.L, seed);

  double *f_wav, *f_scal;
  s2let_allocate_f_wav_real(&f_wav, &f_scal, &parameters);
  s2let_analysis_lm2wav_real(f_wav, f_scal, flm, &parameters);
  s2let_synthesis_wav2lm_real(flm_rec, f_wav, f_scal, &parameters);

  for (int i = 0; i < parameters.L * parameters.L; i += 1) {
    assert_float_equal(creal(flm[i]), creal(flm_rec[i]), 1e-12);
    assert_float_equal(cimag(flm[i]), cimag(flm_rec[i]), 1e-12);
  }

  free(flm);
  free(flm_rec);
  free(f_wav);
  free(f_scal);
}

static s2let_parameters_t state(int spin, int upsample) {
  const s2let_parameters_t parameters = {
      .L = 16,
      .N = 4,
      .B = 2,
      .J_min = 2,
      .spin = spin,
      .original_spin = 0,
      .upsample = upsample,
      .dl_method = SSHT_DL_RISBO,
      .reality = 0,
      .sampling_scheme = S2LET_SAMPLING_MW,
  };
  return parameters;
};

int main(void) {
  const s2let_parameters_t spinzero = state(0, 1);
  const s2let_parameters_t spinzero_noupsample = state(0, 0);
  const s2let_parameters_t spintwo = state(2, 1);
  const s2let_parameters_t spintwo_noupsample = state(2, 0);

#define S2LET_TESTS(name)                                                              \
  cmocka_unit_test_prestate(name, (void *)&spinzero_noupsample),                       \
      cmocka_unit_test_prestate(name, (void *)&spinzero),                              \
      cmocka_unit_test_prestate(name, (void *)&spintwo),                               \
      cmocka_unit_test_prestate(name, (void *)&spintwo_noupsample)

  const struct CMUnitTest tests[] = {
      S2LET_TESTS(test_transform_axisym_lm_wav),
      S2LET_TESTS(test_transform_axisym_lm_wav_multires),
      S2LET_TESTS(test_wav_transform_wavlm_manual),
      S2LET_TESTS(test_wav_transform_harmonic),
      S2LET_TESTS(test_transform_axisym_wav),
      S2LET_TESTS(test_transform_axisym_wav_multires),
      S2LET_TESTS(test_transform_axisym_wav_real),
      S2LET_TESTS(test_transform_axisym_wav_multires_real),
      S2LET_TESTS(test_wav_transform_mw),
      S2LET_TESTS(test_wav_transform_mw_real),
      S2LET_TESTS(test_wav_transform_mwss),
      S2LET_TESTS(test_wav_transform_mwss_real),
      S2LET_TESTS(test_wav_transform_lm2wav),
      S2LET_TESTS(test_wav_transform_lm2wav_real),
  };
#undef S2LET_TESTS

  return cmocka_run_group_tests(tests, NULL, NULL);
}
