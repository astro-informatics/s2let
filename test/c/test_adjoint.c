#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>

#include "s2let/s2let.h"

#include <cmocka.h>

const int seed = 1;

static complex double
dot_complex(complex double *v1, complex double *v2, int N_length) {
  complex long double dot_product = 0.0;
  for (int i = 0; i < N_length; i++)
    dot_product += conj(v1[i]) * v2[i];
  return dot_product;
}

static complex double dot_real(double *v1, double *v2, int N_length) {

  long double dot_product = 0.0;

  for (int i = 0; i < N_length; i++)
    dot_product += v1[i] * v2[i];
  return dot_product;
}

void test_wav_analysis_adjoint_lm_lmn(void **state) {
  s2let_parameters_t *const parameters = (s2let_parameters_t *)*state;

  complex double *flm, *flm_rec;
  s2let_allocate_lm(&flm, parameters->L);
  s2let_allocate_lm(&flm_rec, parameters->L);

  complex double *f_wav_lmn, *f_scal_lm, *f_wav_lmn_rec, *f_scal_lm_rec;
  s2let_allocate_lmn_f_wav(&f_wav_lmn, &f_scal_lm, parameters);
  s2let_allocate_lmn_f_wav(&f_wav_lmn_rec, &f_scal_lm_rec, parameters);

  s2let_lm_random_flm(flm, parameters->L, parameters->spin, seed);
  s2let_lm_random_flm(flm_rec, parameters->L, parameters->spin, seed);

  complex double *wav_lm;
  double *scal_l;
  s2let_tiling_wavelet_allocate(&wav_lm, &scal_l, parameters);
  s2let_tiling_wavelet(wav_lm, scal_l, parameters);

  // Perform wavelet analysis to make f_wac_rec amd f_scal_rec
  s2let_analysis_lm2lmn(
      f_wav_lmn_rec, f_scal_lm_rec, flm_rec, wav_lm, scal_l, parameters);

  // Perform wavelet analysis from scratch with all signals given on the sphere (MW
  s2let_analysis_lm2lmn(f_wav_lmn, f_scal_lm, flm, wav_lm, scal_l, parameters);
  s2let_analysis_adjoint_lmn2lm(
      flm_rec, f_wav_lmn_rec, f_scal_lm_rec, wav_lm, scal_l, parameters);

  // compute dot products of arrays
  complex double error =
      dot_complex(f_wav_lmn_rec, f_wav_lmn, s2let_n_lmn_wav(parameters)) +
      dot_complex(f_scal_lm_rec, f_scal_lm, s2let_n_lm_scal(parameters)) -
      dot_complex(flm_rec, flm, parameters->L * parameters->L);

  assert_float_equal(creal(error), 0, 1e-12);
  assert_float_equal(cimag(error), 0, 1e-12);

  free(flm);
  free(flm_rec);
  free(f_wav_lmn);
  free(f_scal_lm);
  free(f_wav_lmn_rec);
  free(f_scal_lm_rec);
}

void test_wav_synthesis_adjoint_lm_lmn(void **state) {
  s2let_parameters_t *const parameters = (s2let_parameters_t *)*state;

  complex double *flm, *flm_rec;
  s2let_allocate_lm(&flm, parameters->L);
  s2let_allocate_lm(&flm_rec, parameters->L);

  complex double *f_wav_lmn, *f_scal_lm, *f_wav_lmn_rec, *f_scal_lm_rec;
  s2let_allocate_lmn_f_wav(&f_wav_lmn, &f_scal_lm, parameters);
  s2let_allocate_lmn_f_wav(&f_wav_lmn_rec, &f_scal_lm_rec, parameters);

  s2let_lm_random_flm(flm, parameters->L, parameters->spin, seed);
  s2let_lm_random_flm(flm_rec, parameters->L, parameters->spin, seed);

  complex double *wav_lm;
  double *scal_l;
  s2let_tiling_wavelet_allocate(&wav_lm, &scal_l, parameters);
  s2let_tiling_wavelet(wav_lm, scal_l, parameters);

  s2let_analysis_lm2lmn(
      f_wav_lmn_rec, f_scal_lm_rec, flm_rec, wav_lm, scal_l, parameters);

  s2let_synthesis_adjoint_lm2lmn(f_wav_lmn, f_scal_lm, flm, wav_lm, scal_l, parameters);

  s2let_synthesis_lmn2lm(
      flm_rec, f_wav_lmn_rec, f_scal_lm_rec, wav_lm, scal_l, parameters);

  complex double error =
      dot_complex(f_wav_lmn_rec, f_wav_lmn, s2let_n_lmn_wav(parameters)) +
      dot_complex(f_scal_lm_rec, f_scal_lm, s2let_n_lm_scal(parameters)) -
      dot_complex(flm_rec, flm, parameters->L * parameters->L);

  assert_float_equal(creal(error), 0, 1e-12);
  assert_float_equal(cimag(error), 0, 1e-12);

  free(flm);
  free(flm_rec);
  free(f_wav_lmn);
  free(f_scal_lm);
  free(f_wav_lmn_rec);
  free(f_scal_lm_rec);
}

void test_wav_analysis_adjoint_mw(void **state) {
  s2let_parameters_t *const parameters = (s2let_parameters_t *)*state;

  complex double *f, *f_rec, *flm, *flm_rec;
  s2let_allocate_lm(&flm, parameters->L);
  s2let_allocate_lm(&flm_rec, parameters->L);
  s2let_allocate_mw(&f, parameters->L);
  s2let_allocate_mw(&f_rec, parameters->L);

  complex double *f_wav, *f_scal, *f_wav_rec, *f_scal_rec;
  s2let_allocate_f_wav(&f_wav, &f_scal, parameters);
  s2let_allocate_f_wav(&f_wav_rec, &f_scal_rec, parameters);

  s2let_lm_random_flm(flm, parameters->L, parameters->spin, seed);
  s2let_lm_random_flm(flm_rec, parameters->L, parameters->spin, seed);

  ssht_core_mw_inverse_sov_sym(
      f,
      flm,
      parameters->L,
      parameters->spin,
      parameters->dl_method,
      parameters->verbosity);
  ssht_core_mw_inverse_sov_sym(
      f_rec,
      flm_rec,
      parameters->L,
      parameters->spin,
      parameters->dl_method,
      parameters->verbosity);

  s2let_analysis_px2wav(f_wav_rec, f_scal_rec, f_rec, parameters);

  s2let_analysis_px2wav(f_wav, f_scal, f, parameters);

  s2let_analysis_adjoint_wav2px(f_rec, f_wav_rec, f_scal_rec, parameters);

  const complex double error =
      dot_complex(f_wav_rec, f_wav, s2let_n_wav(parameters)) +
      dot_complex(f_scal_rec, f_scal, s2let_n_scal(parameters)) -
      dot_complex(f_rec, f, parameters->L * (2 * parameters->L - 1));

  assert_float_equal(creal(error), 0, 1e-11);
  assert_float_equal(cimag(error), 0, 1e-12);

  free(f);
  free(f_rec);
  free(flm);
  free(flm_rec);
  free(f_wav);
  free(f_scal);
  free(f_wav_rec);
  free(f_scal_rec);
}

void test_wav_analysis_adjoint_mw_real(void **state) {
  s2let_parameters_t parameters = *(s2let_parameters_t *)*state;
  parameters.reality = 1;

  double *f, *f_rec;
  complex double *flm, *flm_rec;
  s2let_allocate_lm(&flm, parameters.L);
  s2let_allocate_lm(&flm_rec, parameters.L);
  s2let_allocate_mw_real(&f, parameters.L);
  s2let_allocate_mw_real(&f_rec, parameters.L);

  double *f_wav, *f_scal, *f_wav_rec, *f_scal_rec;
  s2let_allocate_f_wav_real(&f_wav, &f_scal, &parameters);
  s2let_allocate_f_wav_real(&f_wav_rec, &f_scal_rec, &parameters);

  s2let_lm_random_flm_real(flm, parameters.L, seed);
  s2let_lm_random_flm_real(flm_rec, parameters.L, seed);

  ssht_core_mw_inverse_sov_sym_real(
      f, flm, parameters.L, parameters.dl_method, parameters.verbosity);
  ssht_core_mw_inverse_sov_sym_real(
      f_rec, flm_rec, parameters.L, parameters.dl_method, parameters.verbosity);

  s2let_analysis_px2wav_real(f_wav_rec, f_scal_rec, f_rec, &parameters);
  s2let_analysis_px2wav_real(f_wav, f_scal, f, &parameters);

  s2let_analysis_adjoint_wav2px_real(f_rec, f_wav_rec, f_scal_rec, &parameters);

  const double error = dot_real(f_wav_rec, f_wav, s2let_n_wav(&parameters)) +
                       dot_real(f_scal_rec, f_scal, s2let_n_scal(&parameters)) -
                       dot_real(f_rec, f, parameters.L * (2 * parameters.L - 1));

  assert_float_equal(error, 0, 1e-11);

  free(f);
  free(f_rec);
  free(flm);
  free(flm_rec);
  free(f_wav);
  free(f_scal);
  free(f_wav_rec);
  free(f_scal_rec);
}

void test_wav_synthesis_adjoint_lm2wav(void **state) {
  skip();
  s2let_parameters_t *const parameters = (s2let_parameters_t *)*state;

  complex double *f, *f_rec, *flm, *flm_rec;
  s2let_allocate_lm(&flm, parameters->L);
  s2let_allocate_lm(&flm_rec, parameters->L);
  s2let_allocate_mw(&f, parameters->L);
  s2let_allocate_mw(&f_rec, parameters->L);

  complex double *f_wav, *f_scal, *f_wav_rec, *f_scal_rec;
  s2let_allocate_f_wav(&f_wav, &f_scal, parameters);
  s2let_allocate_f_wav(&f_wav_rec, &f_scal_rec, parameters);

  s2let_lm_random_flm(flm, parameters->L, parameters->spin, seed);
  s2let_lm_random_flm(flm_rec, parameters->L, parameters->spin, seed);

  ssht_core_mw_inverse_sov_sym(
      f,
      flm,
      parameters->L,
      parameters->spin,
      parameters->dl_method,
      parameters->verbosity);
  ssht_core_mw_inverse_sov_sym(
      f_rec,
      flm_rec,
      parameters->L,
      parameters->spin,
      parameters->dl_method,
      parameters->verbosity);

  s2let_analysis_px2wav(f_wav_rec, f_scal_rec, f_rec, parameters);
  s2let_synthesis_adjoint_lm2wav(f_wav, f_scal, flm, parameters);
  s2let_synthesis_wav2lm(flm_rec, f_wav_rec, f_scal_rec, parameters);

  const complex double error =
      dot_complex(f_wav_rec, f_wav, s2let_n_wav(parameters)) +
      dot_complex(f_scal_rec, f_scal, s2let_n_scal(parameters)) -
      dot_complex(f_rec, f, parameters->L * parameters->L);

  assert_float_equal(creal(error), 0, 1e-12);
  assert_float_equal(cimag(error), 0, 1e-12);

  free(f);
  free(f_rec);
  free(flm);
  free(flm_rec);
  free(f_wav);
  free(f_scal);
  free(f_wav_rec);
  free(f_scal_rec);
}

void test_wav_synthesis_adjoint_mw(void **state) {
  s2let_parameters_t *const parameters = (s2let_parameters_t *)*state;

  complex double *f, *f_rec, *flm, *flm_rec;
  s2let_allocate_lm(&flm, parameters->L);
  s2let_allocate_lm(&flm_rec, parameters->L);
  s2let_allocate_mw(&f, parameters->L);
  s2let_allocate_mw(&f_rec, parameters->L);

  complex double *f_wav, *f_scal, *f_wav_rec, *f_scal_rec;
  s2let_allocate_f_wav(&f_wav, &f_scal, parameters);
  s2let_allocate_f_wav(&f_wav_rec, &f_scal_rec, parameters);

  s2let_lm_random_flm(flm, parameters->L, parameters->spin, seed);
  s2let_lm_random_flm(flm_rec, parameters->L, parameters->spin, seed);

  ssht_core_mw_inverse_sov_sym(
      f,
      flm,
      parameters->L,
      parameters->spin,
      parameters->dl_method,
      parameters->verbosity);
  ssht_core_mw_inverse_sov_sym(
      f_rec,
      flm_rec,
      parameters->L,
      parameters->spin,
      parameters->dl_method,
      parameters->verbosity);

  s2let_analysis_lm2wav(f_wav_rec, f_scal_rec, flm_rec, parameters);

  s2let_synthesis_adjoint_px2wav(f_wav, f_scal, f, parameters);
  s2let_synthesis_wav2px(f_rec, f_wav_rec, f_scal_rec, parameters);

  const complex double error =
      dot_complex(f_wav_rec, f_wav, s2let_n_wav(parameters)) +
      dot_complex(f_scal_rec, f_scal, s2let_n_scal(parameters)) -
      dot_complex(f_rec, f, parameters->L * (2 * parameters->L - 1));

  assert_float_equal(creal(error), 0, 1e-12);
  assert_float_equal(cimag(error), 0, 1e-12);

  free(f);
  free(f_rec);
  free(flm);
  free(flm_rec);
  free(f_wav);
  free(f_scal);
  free(f_wav_rec);
  free(f_scal_rec);
}

void test_wav_synthesis_adjoint_mw_real(void **state) {
  s2let_parameters_t parameters = *(s2let_parameters_t *)*state;
  parameters.reality = 1;

  double *f, *f_rec;
  complex double *flm, *flm_rec;
  s2let_allocate_lm(&flm, parameters.L);
  s2let_allocate_lm(&flm_rec, parameters.L);
  s2let_allocate_mw_real(&f, parameters.L);
  s2let_allocate_mw_real(&f_rec, parameters.L);

  double *f_wav, *f_scal, *f_wav_rec, *f_scal_rec;
  s2let_allocate_f_wav_real(&f_wav, &f_scal, &parameters);
  s2let_allocate_f_wav_real(&f_wav_rec, &f_scal_rec, &parameters);

  s2let_lm_random_flm_real(flm, parameters.L, seed);
  s2let_lm_random_flm_real(flm_rec, parameters.L, seed);

  ssht_core_mw_inverse_sov_sym_real(
      f, flm, parameters.L, parameters.dl_method, parameters.verbosity);
  ssht_core_mw_inverse_sov_sym_real(
      f_rec, flm_rec, parameters.L, parameters.dl_method, parameters.verbosity);

  s2let_analysis_px2wav_real(f_wav_rec, f_scal_rec, f_rec, &parameters);

  s2let_synthesis_adjoint_px2wav_real(f_wav, f_scal, f, &parameters);

  s2let_synthesis_wav2px_real(f_rec, f_wav_rec, f_scal_rec, &parameters);

  const double error = dot_real(f_wav_rec, f_wav, s2let_n_wav(&parameters)) +
                       dot_real(f_scal_rec, f_scal, s2let_n_scal(&parameters)) -
                       dot_real(f_rec, f, parameters.L * (2 * parameters.L - 1));

  assert_float_equal(error, 0, 1e-11);

  free(f);
  free(f_rec);
  free(flm);
  free(flm_rec);
  free(f_wav);
  free(f_scal);
  free(f_wav_rec);
  free(f_scal_rec);
}

void gen_flmn_complex(
    complex double *flmn, const so3_parameters_t *parameters, int seed) {
  int L0, L, N;
  int i, el, m, n, n_start, n_stop, n_inc, ind;

  L0 = parameters->L0;
  L = parameters->L;
  N = parameters->N;

  for (i = 0; i < (2 * N - 1) * L * L; ++i)
    flmn[i] = 0.0;

  switch (parameters->n_mode) {
  case SO3_N_MODE_ALL:
    n_start = -N + 1;
    n_stop = N - 1;
    n_inc = 1;
    break;
  case SO3_N_MODE_EVEN:
    n_start = ((N - 1) % 2 == 0) ? -N + 1 : -N + 2;
    n_stop = ((N - 1) % 2 == 0) ? N - 1 : N - 2;
    n_inc = 2;
    break;
  case SO3_N_MODE_ODD:
    n_start = ((N - 1) % 2 != 0) ? -N + 1 : -N + 2;
    n_stop = ((N - 1) % 2 != 0) ? N - 1 : N - 2;
    n_inc = 2;
    break;
  case SO3_N_MODE_MAXIMUM:
    n_start = -N + 1;
    n_stop = N - 1;
    n_inc = 2 * N - 2;
    break;
  default:
    SO3_ERROR_GENERIC("Invalid n-mode.");
  }

  for (n = n_start; n <= n_stop; n += n_inc) {
    for (el = MAX(L0, abs(n)); el < L; ++el) {
      for (m = -el; m <= el; ++m) {
        so3_sampling_elmn2ind(&ind, el, m, n, parameters);
        flmn[ind] = (2.0 * ran2_dp(seed) - 1.0) + I * (2.0 * ran2_dp(seed) - 1.0);
      }
    }
  }
}

void test_wav_so3_forward_adjoint(void **state) {
  const s2let_parameters_t *parameters = (const s2let_parameters_t *)*state;

  so3_parameters_t so3_parameters;
  fill_so3_parameters(&so3_parameters, parameters);

  complex double *f, *f_rec, *flmn, *flmn_rec;
  flmn = malloc((2 * parameters->N - 1) * parameters->L * parameters->L * sizeof *flmn);
  SO3_ERROR_MEM_ALLOC_CHECK(flmn);
  flmn_rec = malloc(
      (2 * parameters->N - 1) * parameters->L * parameters->L * sizeof *flmn_rec);
  SO3_ERROR_MEM_ALLOC_CHECK(flmn_rec);

  f = malloc(
      (2 * parameters->L - 1) * (parameters->L) * (2 * parameters->N - 1) * sizeof *f);
  SO3_ERROR_MEM_ALLOC_CHECK(f);
  f_rec = malloc(
      (2 * parameters->L - 1) * (parameters->L) * (2 * parameters->N - 1) *
      sizeof *f_rec);
  SO3_ERROR_MEM_ALLOC_CHECK(f_rec);

  gen_flmn_complex(flmn, &so3_parameters, seed);
  gen_flmn_complex(flmn_rec, &so3_parameters, seed);

  so3_core_inverse_direct(f_rec, flmn_rec, &so3_parameters);

  so3_adjoint_forward_direct(f, flmn, &so3_parameters);
  so3_core_forward_direct(flmn_rec, f_rec, &so3_parameters);

  const complex double error =
      dot_complex(
          f_rec,
          f,
          (2 * parameters->L - 1) * (parameters->L) * (2 * parameters->N - 1)) -
      dot_complex(
          flmn_rec, flmn, parameters->L * parameters->L * (2 * parameters->N - 1));

  assert_float_equal(creal(error), 0, 1e-12);
  assert_float_equal(cimag(error), 0, 1e-12);

  free(f);
  free(f_rec);
  free(flmn);
  free(flmn_rec);
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
  const s2let_parameters_t spintwo = state(2, 1);
  const s2let_parameters_t noupsample = state(2, 0);
  const struct CMUnitTest tests[] = {
      cmocka_unit_test_prestate(test_wav_analysis_adjoint_lm_lmn, (void *)&spinzero),
      cmocka_unit_test_prestate(test_wav_analysis_adjoint_lm_lmn, (void *)&spintwo),
      cmocka_unit_test_prestate(test_wav_analysis_adjoint_lm_lmn, (void *)&noupsample),

      cmocka_unit_test_prestate(test_wav_analysis_adjoint_mw, (void *)&spinzero),
      cmocka_unit_test_prestate(test_wav_analysis_adjoint_mw, (void *)&spintwo),
      cmocka_unit_test_prestate(test_wav_analysis_adjoint_mw, (void *)&noupsample),

      cmocka_unit_test_prestate(test_wav_analysis_adjoint_mw_real, (void *)&spinzero),
      cmocka_unit_test_prestate(test_wav_analysis_adjoint_mw_real, (void *)&spintwo),
      cmocka_unit_test_prestate(test_wav_analysis_adjoint_mw_real, (void *)&noupsample),

      cmocka_unit_test_prestate(test_wav_synthesis_adjoint_lm_lmn, (void *)&spinzero),
      cmocka_unit_test_prestate(test_wav_synthesis_adjoint_lm_lmn, (void *)&spintwo),
      cmocka_unit_test_prestate(test_wav_synthesis_adjoint_lm_lmn, (void *)&noupsample),

      cmocka_unit_test_prestate(test_wav_synthesis_adjoint_lm2wav, (void *)&spinzero),
      cmocka_unit_test_prestate(test_wav_synthesis_adjoint_lm2wav, (void *)&spintwo),
      cmocka_unit_test_prestate(test_wav_synthesis_adjoint_lm2wav, (void *)&noupsample),

      cmocka_unit_test_prestate(test_wav_synthesis_adjoint_mw, (void *)&spinzero),
      cmocka_unit_test_prestate(test_wav_synthesis_adjoint_mw, (void *)&spintwo),
      cmocka_unit_test_prestate(test_wav_synthesis_adjoint_mw, (void *)&noupsample),

      cmocka_unit_test_prestate(test_wav_synthesis_adjoint_mw_real, (void *)&spinzero),
      cmocka_unit_test_prestate(test_wav_synthesis_adjoint_mw_real, (void *)&spintwo),
      cmocka_unit_test_prestate(
          test_wav_synthesis_adjoint_mw_real, (void *)&noupsample),

      cmocka_unit_test_prestate(test_wav_so3_forward_adjoint, (void *)&spinzero),
      cmocka_unit_test_prestate(test_wav_so3_forward_adjoint, (void *)&spintwo),
      cmocka_unit_test_prestate(test_wav_so3_forward_adjoint, (void *)&noupsample),
  };

  return cmocka_run_group_tests(tests, NULL, NULL);
}
