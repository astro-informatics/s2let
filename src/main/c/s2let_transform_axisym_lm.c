// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include <complex.h>
#include <math.h>
#include <stdlib.h>

#include "s2let/s2let.h"

static inline int lm2ind(int el, int m) { return el * el + el + m; }

/*!
 * Allocates arrays for wavelet transform in harmonic space.
 *
 * \param[out]  f_wav_lm Harmonic coefficients of the wavelet contributions.
 *                       Each wavelet has size L*L and there are (J-J_min+1)
 *                       scales.
 * \param[out]  f_scal_lm Harmonic coefficients of the scaling contribution (L*L).
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_lm_allocate_f_wav(
    complex double **f_wav_lm,
    complex double **f_scal_lm,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int J = s2let_j_max(parameters);
  *f_wav_lm = calloc((J + 1 - J_min) * L * L, sizeof **f_wav_lm);
  *f_scal_lm = calloc(L * L, sizeof **f_scal_lm);
}

/*!
 * Allocates multiresolution arrays for wavelet transforms in harmonic space.
 *
 * \param[out]  f_wav_lm Harmonic coefficients of the wavelet contributions. The size of
 * each wavelet depends on its band-limit. \param[out]  f_scal_lm Harmonic coefficients
 * of the scaling contribution. \param[in]  B Wavelet parameter. \param[in]  L Angular
 * harmonic band-limit. \param[in]  J_min First wavelet scale to be used. \retval none
 */
void s2let_transform_axisym_lm_allocate_f_wav_multires(
    complex double **f_wav_lm,
    complex double **f_scal_lm,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int J = s2let_j_max(parameters);
  int j, bandlimit, total = 0;
  for (j = J_min; j <= J; ++j) {
    bandlimit = MIN(s2let_bandlimit(j, parameters), L);
    total += bandlimit * bandlimit;
  }
  *f_wav_lm = calloc(total, sizeof **f_wav_lm);
  bandlimit = MIN(s2let_bandlimit(J_min - 1, parameters), L);
  *f_scal_lm = calloc(bandlimit * bandlimit, sizeof **f_scal_lm);
}

/*!
 * Allocates arrays for the kernels of the wavelets and the scaling functions.
 *
 * \param[out]  wav_lm Wavelet kernels.
 * \param[out]  scal_lm Scaling function kernels.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_transform_axisym_lm_allocate_wav(
    double **wav_lm, double **scal_lm, const s2let_parameters_t *parameters) {
  int L = parameters->L;

  int J = s2let_j_max(parameters);
  *wav_lm = (double *)calloc((J + 1) * L, sizeof(double));
  *scal_lm = (double *)calloc(L, sizeof(double));
}

/*!
 * Computes the kernels of the wavelets and the scaling functions.
 *
 * \param[out]  wav_lm Wavelet kernels.
 * \param[out]  scal_lm Scaling function kernels.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_lm_wav(
    double *wav_lm, double *scal_lm, const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int j, l;
  int J = s2let_j_max(parameters);
  // int J_min = 0;
  // int l_min = s2let_transform_axisym_el_min(B, J_min);
  double k0;
  double *kappa, *kappa0;
  s2let_tiling_axisym_allocate(&kappa, &kappa0, parameters);
  s2let_tiling_axisym(kappa, kappa0, parameters);

  for (j = J_min; j <= J; j++) {
    for (l = 0; l < L; l++) {
      k0 = sqrt((2 * l + 1) / (4.0 * PI)) * kappa[l + j * L];
      wav_lm[j * L + l] = k0;
    }
  }
  for (l = 0; l < L; l++) {
    k0 = sqrt((2 * l + 1) / (4.0 * PI)) * kappa0[l];
    scal_lm[l] = k0;
  }

  free(kappa);
  free(kappa0);
}

/*!
 * Spherical wavelets : full resolution analysis in harmonic space.
 * Perform wavelet transform from precomputed kernels and gives the harmonic
 * coefficients.
 *
 * \param[out]  f_wav_lm Wavelet transform (harmonic coefficients of wavelet
 * contribution). \param[out]  f_scal_lm Wavelet transform (harmonic coefficients of
 * scaling contribution). \param[in]  flm Spherical harmonic decomposition of input
 * function. \param[in]  wav_lm Wavelet kernels. \param[in]  scal_lm Scaling function
 * kernels. \param[in]  B Wavelet parameter. \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_lm_wav_analysis(
    complex double *f_wav_lm,
    complex double *f_scal_lm,
    const complex double *flm,
    const double *wav_lm,
    const double *scal_lm,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int offset, j, l, m;
  int J = s2let_j_max(parameters);
  double wav0, scal0;
  // int l_min = s2let_transform_axisym_el_min(B, J_min);

  offset = 0;
  for (j = J_min; j <= J; j++) {
    for (l = 0; l < L; l++) {
      wav0 = sqrt((4.0 * PI) / (2.0 * l + 1.0)) * wav_lm[j * L + l];
      for (m = -l; m <= l; m++) {
        f_wav_lm[offset + l * l + l + m] = flm[lm2ind(l, m)] * wav0;
      }
    }
    offset += L * L;
  }
  for (l = 0; l < L; l++) {
    scal0 = sqrt((4.0 * PI) / (2.0 * l + 1.0)) * scal_lm[l];
    for (m = -l; m <= l; m++) {
      f_scal_lm[lm2ind(l, m)] = flm[lm2ind(l, m)] * scal0;
    }
  }
}

/*!
 * Spherical wavelets : full resolution synthesis in harmonic space.
 * Perform wavelet transform in harmonic space from precomputed kernels and gives
 * harmonic coefficients.
 *
 * \param[out]  flm Spherical harmonic decomposition of input function.
 * \param[in]  f_wav_lm Wavelet transform (harmonic coefficients of wavelet
 * contribution). \param[in]  f_scal_lm Wavelet transform (harmonic coefficients of
 * scaling contribution). \param[in]  wav_lm Wavelet kernels. \param[in]  scal_lm
 * Scaling function kernels. \param[in]  B Wavelet parameter. \param[in]  L Angular
 * harmonic band-limit. \param[in]  J_min First wavelet scale to be used. \retval none
 */
void s2let_transform_axisym_lm_wav_synthesis(
    complex double *flm,
    const complex double *f_wav_lm,
    const complex double *f_scal_lm,
    const double *wav_lm,
    const double *scal_lm,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int offset, j, l, m;
  int J = s2let_j_max(parameters);
  double wav0, scal0;
  // int l_min = s2let_transform_axisym_el_min(B, J_min);

  offset = 0;
  for (j = J_min; j <= J; j++) {
    for (l = 0; l < L; l++) {
      wav0 = sqrt((4.0 * PI) / (2.0 * l + 1.0)) * wav_lm[j * L + l];
      for (m = -l; m <= l; m++) {
        flm[lm2ind(l, m)] += f_wav_lm[offset + l * l + l + m] * wav0;
      }
    }
    offset += L * L;
  }
  for (l = 0; l < L; l++) {
    scal0 = sqrt((4.0 * PI) / (2.0 * l + 1.0)) * scal_lm[l];
    for (m = -l; m <= l; m++) {
      flm[lm2ind(l, m)] += f_scal_lm[lm2ind(l, m)] * scal0;
    }
  }
}

/*!
 * Spherical wavelets : multiresolution analysis in harmonic space.
 * Perform multiresolution wavelet transform in harmonic space from precomputed kernels
 * and gives harmonic coefficients.
 *
 * \param[out]  f_wav_lm Wavelet transform (SHA of wavelet contribution).
 * \param[out]  f_scal_lm Wavelet transform (SHA of scaling contribution).
 * \param[in]  flm Spherical harmonic decomposition of input function.
 * \param[in]  wav_lm Wavelet kernels.
 * \param[in]  scal_lm Scaling function kernels.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_lm_wav_analysis_multires(
    complex double *f_wav_lm,
    complex double *f_scal_lm,
    const complex double *flm,
    const double *wav_lm,
    const double *scal_lm,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int bandlimit, offset, j, l, m;
  int J = s2let_j_max(parameters);
  double wav0, scal0;

  offset = 0;
  for (j = J_min; j <= J; j++) {
    bandlimit = MIN(s2let_bandlimit(j, parameters), L);
    for (l = 0; l < bandlimit; l++) {
      wav0 = sqrt((4.0 * PI) / (2.0 * l + 1.0)) * wav_lm[j * L + l];
      for (m = -l; m <= l; m++) {
        f_wav_lm[offset + l * l + l + m] = flm[lm2ind(l, m)] * wav0;
      }
    }
    offset += bandlimit * bandlimit;
  }
  bandlimit = MIN(s2let_bandlimit(J_min - 1, parameters), L);
  for (l = 0; l < bandlimit; l++) {
    scal0 = sqrt((4.0 * PI) / (2.0 * l + 1.0)) * scal_lm[l];
    for (m = -l; m <= l; m++) {
      f_scal_lm[lm2ind(l, m)] = flm[lm2ind(l, m)] * scal0;
    }
  }
}

/*!
 * Spherical wavelets : multiresolution synthesis in harmonic space.
 * Perform multiresolution wavelet transform in harmonic space from precomputed kernels
 * and gives harmonic coefficients.
 *
 * \param[out]  flm Spherical harmonic decomposition of input function.
 * \param[in]  f_wav_lm Wavelet transform (SHA of wavelet contribution).
 * \param[in]  f_scal_lm Wavelet transform (SHA of scaling contribution).
 * \param[in]  wav_lm Wavelet kernels.
 * \param[in]  scal_lm Scaling function kernels.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_lm_wav_synthesis_multires(
    complex double *flm,
    const complex double *f_wav_lm,
    const complex double *f_scal_lm,
    const double *wav_lm,
    const double *scal_lm,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int bandlimit, offset, j, l, m;
  int J = s2let_j_max(parameters);
  double wav0, scal0;
  // int l_min = s2let_transform_axisym_el_min(B, J_min);

  offset = 0;
  for (j = J_min; j <= J; j++) {
    bandlimit = MIN(s2let_bandlimit(j, parameters), L);
    for (l = 0; l < bandlimit; l++) {
      wav0 = sqrt((4.0 * PI) / (2.0 * l + 1.0)) * wav_lm[j * L + l];
      for (m = -l; m <= l; m++) {
        flm[lm2ind(l, m)] += f_wav_lm[offset + l * l + l + m] * wav0;
      }
    }
    offset += bandlimit * bandlimit;
  }
  bandlimit = MIN(s2let_bandlimit(J_min - 1, parameters), L);
  for (l = 0; l < bandlimit; l++) {
    scal0 = sqrt((4.0 * PI) / (2.0 * l + 1.0)) * scal_lm[l];
    for (m = -l; m <= l; m++) {
      flm[lm2ind(l, m)] += f_scal_lm[lm2ind(l, m)] * scal0;
    }
  }
}
