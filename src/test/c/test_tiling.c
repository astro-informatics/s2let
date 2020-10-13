#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>

#include "s2let.h"

#include <cmocka.h>

const int seed = 1;

/*! Identity relation of the wavelet tiling in harmonic space */
void test_tiling_axisym(void **state) {
  const s2let_parameters_t *const parameters = (const s2let_parameters_t *)*state;

  double *kappa, *kappa0;

  // Allocate the kernels corresponding to the parameters B, L
  s2let_tiling_axisym_allocate(&kappa, &kappa0, parameters);

  // Construct the tiling of harmonic space
  s2let_tiling_axisym(kappa, kappa0, parameters);

  // Check that they recover the identity relation,
  // ensuring exactness of the wavelet transform.
  double res = s2let_tiling_axisym_check_identity(kappa, kappa0, parameters);
  assert_float_equal(res, 0, 1e-12);

  free(kappa);
  free(kappa0);
}

/*! Identity relation of the directionality components for tiling in harmonic space. */
void test_tiling_direction(void **state) {
  const s2let_parameters_t *const parameters = (const s2let_parameters_t *)*state;

  complex double *s_elm;
  double error;

  // Allocate space for the harmonic coefficients
  s2let_tiling_direction_allocate(&s_elm, parameters);

  // Construct the harmonic coefficients
  s2let_tiling_direction(s_elm, parameters);

  // Check that they recover the identity relation,
  // ensuring exactness of the wavelet transform.
  error = s2let_tiling_direction_check_identity(s_elm, parameters);
  assert_float_equal(error, 0, 1e-12);

  free(s_elm);
}

/*! Identity relation of the directional wavelets for tiling in harmonic space. */
void test_tiling_wavelet(void **state) {
  const s2let_parameters_t *const parameters = (const s2let_parameters_t *)*state;

  complex double *phi;
  double *psi;
  double error;

  // Allocate space for the harmonic coefficients
  s2let_tiling_wavelet_allocate(&phi, &psi, parameters);

  // Construct the harmonic coefficients
  s2let_tiling_wavelet(phi, psi, parameters);

  // Check that they recover the identity relation,
  // ensuring exactness of the wavelet transform.
  error = s2let_tiling_wavelet_check_identity(phi, psi, parameters);
  assert_float_equal(error, 0, 1e-12);

  free(phi);
  free(psi);
}

void test_binomial_coefficient(void **state) {
  for (int n = 1; n <= 62; ++n)
    for (int k = 0; k <= n / 2; ++k)
      assert_float_equal(
          binomial_coefficient(n, k, 0), binomial_coefficient(n, k, 1), 1e-12);
}

int main(void) {
  const s2let_parameters_t state = {
      .L = 16,
      .N = 4,
      .B = 2,
      .J_min = 2,
      .spin = 0,
      .original_spin = 0,
  };

  const struct CMUnitTest tests[] = {
      cmocka_unit_test_prestate(test_tiling_axisym, (void *)&state),
      cmocka_unit_test_prestate(test_tiling_direction, (void *)&state),
      cmocka_unit_test_prestate(test_tiling_wavelet, (void *)&state),
      cmocka_unit_test(test_binomial_coefficient),
  };

  return cmocka_run_group_tests(tests, NULL, NULL);
}
