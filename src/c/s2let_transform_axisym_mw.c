// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include <complex.h>
#include <math.h>
#include <ssht/ssht.h>
#include <stdlib.h>

#include "s2let/s2let.h"

/*!
 * Allocates arrays for final wavelets and scaling functions in pixel space (MW
 * sampling).
 *
 * \param[out]  f_wav Array of wavelets maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_allocate_mw_f_wav(
    complex double **f_wav,
    complex double **f_scal,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int J = s2let_j_max(parameters);
  *f_wav = calloc((J + 1 - J_min) * L * (2 * L - 1), sizeof **f_wav);
  *f_scal = calloc(L * (2 * L - 1), sizeof **f_scal);
}

/*!
 * Allocates arrays for wavelets and scaling functions in pixel space (MW sampling).
 *
 * \param[out]  f_wav Array of wavelets maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_allocate_mw_f_wav_multires(
    complex double **f_wav,
    complex double **f_scal,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int J = s2let_j_max(parameters);
  int j, bandlimit, total = 0;
  for (j = J_min; j <= J; ++j) {
    bandlimit = MIN(s2let_bandlimit(j, parameters), L);
    total += bandlimit * (2 * bandlimit - 1);
  }
  *f_wav = calloc(total, sizeof **f_wav);
  bandlimit = MIN(s2let_bandlimit(J_min - 1, parameters), L);
  *f_scal = calloc(bandlimit * (2 * bandlimit - 1), sizeof **f_scal);
}

/*!
 * Allocates arrays for final wavelets and scaling functions in pixel space (MW
 * sampling).
 *
 * \param[out]  f_wav Array of wavelets maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_allocate_mw_f_wav_real(
    double **f_wav, double **f_scal, const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int J = s2let_j_max(parameters);
  *f_wav = calloc((J + 1 - J_min) * L * (2 * L - 1), sizeof **f_wav);
  *f_scal = calloc(L * (2 * L - 1), sizeof **f_scal);
}

/*!
 * Allocates arrays for multiresolution wavelets and scaling functions in pixel space
 * (MW sampling).
 *
 * \param[out]  f_wav Array of wavelets maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_allocate_mw_f_wav_multires_real(
    double **f_wav, double **f_scal, const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int J = s2let_j_max(parameters);
  int j, bandlimit, total = 0;
  for (j = J_min; j <= J; ++j) {
    bandlimit = MIN(s2let_bandlimit(j, parameters), L);
    total += bandlimit * (2 * bandlimit - 1);
  }
  *f_wav = calloc(total, sizeof **f_wav);
  bandlimit = MIN(s2let_bandlimit(J_min - 1, parameters), L);
  *f_scal = calloc(bandlimit * (2 * bandlimit - 1), sizeof **f_scal);
}

/*!
 * Spherical wavelets : full resolution analysis in real space, MW sampling.
 * Perform wavelet transform in real space (from scratch, gives wavelet maps).
 * Sampling scheme : MW sampling.
 *
 * \param[out]  f_wav Array of wavelets maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling..
 * \param[in]  f Input function (MW sampling)
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_wav_analysis_mw(
    complex double *f_wav,
    complex double *f_scal,
    const complex double *f,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int spin = 0;
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;

  int j, offset, offset_lm;
  int J = s2let_j_max(parameters);
  // int l_min = s2let_transform_axisym_el_min(B, J_min);

  double *wav_lm, *scal_lm;
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, parameters);
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, parameters);

  complex double *flm, *f_wav_lm, *f_scal_lm;
  flm = (complex double *)calloc(L * L, sizeof(complex double));
  s2let_transform_axisym_lm_allocate_f_wav(&f_wav_lm, &f_scal_lm, parameters);

  ssht_core_mw_forward_sov_conv_sym(flm, f, L, spin, dl_method, verbosity);

  s2let_transform_axisym_lm_wav_analysis(
      f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, parameters);

  ssht_core_mw_inverse_sov_sym(f_scal, f_scal_lm, L, spin, dl_method, verbosity);
  offset = 0;
  offset_lm = 0;
  for (j = J_min; j <= J; ++j) {
    int i;

    ssht_core_mw_inverse_sov_sym(
        f_wav + offset, f_wav_lm + offset_lm, L, spin, dl_method, verbosity);

    for (i = 0; i < L * (2 * L - 1); ++i) {
      f_wav[offset + i] /= sqrt(2 * PI);
    }

    offset_lm += L * L;
    offset += L * (2 * L - 1);
  }

  free(flm);
  free(f_scal_lm);
  free(f_wav_lm);
}

/*!
 * Perform wavelet transform adjoint in real space (from scratch, gives wavelet maps).
 * Sampling scheme : MW sampling.
 * Spherical wavelets : analysis adjoint in real space, MW sampling.
 *
 * \param[out]  f Input function (MW sampling)
 * \param[in]  f_wav Array of wavelets maps, MW sampling.
 * \param[in]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_wav_analysis_adjoint_mw(
    complex double *f,
    const complex double *f_wav,
    const complex double *f_scal,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int spin = 0;
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;

  int j, offset, offset_lm;
  int J = s2let_j_max(parameters);
  // int l_min = s2let_transform_axisym_el_min(B, J_min);

  double *wav_lm, *scal_lm;
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, parameters);
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, parameters);

  complex double *flm, *f_wav_norm, *f_wav_lm, *f_scal_lm;
  flm = (complex double *)calloc(L * L, sizeof(complex double));
  s2let_transform_axisym_lm_allocate_f_wav(&f_wav_lm, &f_scal_lm, parameters);

  ssht_adjoint_mw_inverse_sov_sym(f_scal_lm, f_scal, L, spin, dl_method, verbosity);
  offset = 0;
  offset_lm = 0;

  f_wav_norm = malloc(L * (2 * L - 1) * sizeof *f_wav_norm);
  for (j = J_min; j <= J; j++) {
    int i;

    for (i = 0; i < L * (2 * L - 1); ++i) {
      f_wav_norm[i] = f_wav[offset + i] / sqrt(2 * PI);
    }

    ssht_adjoint_mw_inverse_sov_sym(
        f_wav_lm + offset_lm, f_wav_norm, L, spin, dl_method, verbosity);
    offset_lm += L * L;
    offset += L * (2 * L - 1);
  }
  free(f_wav_norm);

  s2let_transform_axisym_lm_wav_synthesis(
      flm, f_wav_lm, f_scal_lm, wav_lm, scal_lm, parameters);

  ssht_adjoint_mw_forward_sov_sym(f, flm, L, spin, dl_method, verbosity);

  free(flm);
  free(f_scal_lm);
  free(f_wav_lm);
}

/*!
 * Perform wavelet transform in real space (from scratch, gives wavelet maps).
 * Sampling scheme : MW sampling.
 * Spherical wavelets : synthesis in real space, MW sampling.
 *
 * \param[out]  f Input function (MW sampling)
 * \param[in]  f_wav Array of wavelets maps, MW sampling.
 * \param[in]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_wav_synthesis_mw(
    complex double *f,
    const complex double *f_wav,
    const complex double *f_scal,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int spin = 0;
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;

  int j, offset, offset_lm;
  int J = s2let_j_max(parameters);
  // int l_min = s2let_transform_axisym_el_min(B, J_min);

  double *wav_lm, *scal_lm;
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, parameters);
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, parameters);

  complex double *flm, *f_wav_norm, *f_wav_lm, *f_scal_lm;
  flm = (complex double *)calloc(L * L, sizeof(complex double));
  s2let_transform_axisym_lm_allocate_f_wav(&f_wav_lm, &f_scal_lm, parameters);

  ssht_core_mw_forward_sov_conv_sym(f_scal_lm, f_scal, L, spin, dl_method, verbosity);
  offset = 0;
  offset_lm = 0;

  f_wav_norm = malloc(L * (2 * L - 1) * sizeof *f_wav_norm);
  for (j = J_min; j <= J; j++) {
    int i;

    for (i = 0; i < L * (2 * L - 1); ++i) {
      f_wav_norm[i] = f_wav[offset + i] * sqrt(2 * PI);
    }

    ssht_core_mw_forward_sov_conv_sym(
        f_wav_lm + offset_lm, f_wav_norm, L, spin, dl_method, verbosity);
    offset_lm += L * L;
    offset += L * (2 * L - 1);
  }
  free(f_wav_norm);

  s2let_transform_axisym_lm_wav_synthesis(
      flm, f_wav_lm, f_scal_lm, wav_lm, scal_lm, parameters);

  ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

  free(flm);
  free(f_scal_lm);
  free(f_wav_lm);
}

/*!
 * Spherical wavelets : full resolution synthesis in real space adjoint, MW sampling.
 * Perform wavelet synthesis adjoint in real space (from scratch, gives wavelet maps).
 * Sampling scheme : MW sampling.
 *
 * \param[out]  f_wav Array of wavelets maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling..
 * \param[in]  f Input function (MW sampling)
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_wav_synthesis_adjoint_mw(
    complex double *f_wav,
    complex double *f_scal,
    const complex double *f,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int spin = 0;
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;

  int j, offset, offset_lm;
  int J = s2let_j_max(parameters);
  // int l_min = s2let_transform_axisym_el_min(B, J_min);

  double *wav_lm, *scal_lm;
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, parameters);
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, parameters);

  complex double *flm, *f_wav_lm, *f_scal_lm;
  flm = (complex double *)calloc(L * L, sizeof(complex double));
  s2let_transform_axisym_lm_allocate_f_wav(&f_wav_lm, &f_scal_lm, parameters);

  ssht_adjoint_mw_inverse_sov_sym(flm, f, L, spin, dl_method, verbosity);

  s2let_transform_axisym_lm_wav_analysis(
      f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, parameters);

  ssht_adjoint_mw_forward_sov_sym(f_scal, f_scal_lm, L, spin, dl_method, verbosity);
  offset = 0;
  offset_lm = 0;
  for (j = J_min; j <= J; ++j) {
    int i;

    ssht_adjoint_mw_forward_sov_sym(
        f_wav + offset, f_wav_lm + offset_lm, L, spin, dl_method, verbosity);

    for (i = 0; i < L * (2 * L - 1); ++i) {
      f_wav[offset + i] *= sqrt(2 * PI);
    }

    offset_lm += L * L;
    offset += L * (2 * L - 1);
  }

  free(flm);
  free(f_scal_lm);
  free(f_wav_lm);
}

/*!
 * Spherical wavelets : multiresolution analysis in real space, MW sampling.
 * Perform multiresolution wavelet transform in real space (from scratch, gives wavelet
 * maps). Sampling scheme : MW sampling.
 *
 * \param[out]  f_wav Array of wavelets maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  f Input function (MW sampling)
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_wav_analysis_mw_multires(
    complex double *f_wav,
    complex double *f_scal,
    const complex double *f,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int spin = 0;
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;

  int bandlimit, j, offset, offset_lm;
  int J = s2let_j_max(parameters);
  // int l_min = s2let_transform_axisym_el_min(B, J_min);

  double *wav_lm, *scal_lm;
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, parameters);
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, parameters);

  complex double *flm, *f_wav_lm, *f_scal_lm;
  flm = (complex double *)calloc(L * L, sizeof(complex double));
  s2let_transform_axisym_lm_allocate_f_wav_multires(&f_wav_lm, &f_scal_lm, parameters);

  ssht_core_mw_forward_sov_conv_sym(flm, f, L, spin, dl_method, verbosity);

  s2let_transform_axisym_lm_wav_analysis_multires(
      f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, parameters);

  bandlimit = MIN(s2let_bandlimit(J_min - 1, parameters), L);
  ssht_core_mw_inverse_sov_sym(
      f_scal, f_scal_lm, bandlimit, spin, dl_method, verbosity);
  offset = 0;
  offset_lm = 0;
  for (j = J_min; j <= J; j++) {
    int i;

    bandlimit = MIN(s2let_bandlimit(j, parameters), L);
    ssht_core_mw_inverse_sov_sym(
        f_wav + offset, f_wav_lm + offset_lm, bandlimit, spin, dl_method, verbosity);

    for (i = 0; i < bandlimit * (2 * bandlimit - 1); ++i) {
      f_wav[offset + i] /= sqrt(2 * PI);
    }

    offset_lm += bandlimit * bandlimit;
    offset += bandlimit * (2 * bandlimit - 1);
  }

  free(flm);
  free(f_scal_lm);
  free(f_wav_lm);
}

/*!
 * Spherical wavelets : multiresolution analysis adjoint in real space, MW sampling.
 * Perform multiresolution wavelet transform adjoint in real space (from scratch, gives
 * wavelet maps). Sampling scheme : MW sampling.
 *
 * \param[out]  f Input function (MW sampling)
 * \param[in]  f_wav Array of wavelets maps, MW sampling.
 * \param[in]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_wav_analysis_adjoint_mw_multires(
    complex double *f,
    const complex double *f_wav,
    const complex double *f_scal,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int spin = 0;
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;

  int bandlimit, j, offset, offset_lm;
  int J = s2let_j_max(parameters);
  // int l_min = s2let_transform_axisym_el_min(B, J_min);

  double *wav_lm, *scal_lm;
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, parameters);
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, parameters);

  complex double *flm, *f_wav_norm, *f_wav_lm, *f_scal_lm;
  flm = (complex double *)calloc(L * L, sizeof(complex double));
  s2let_transform_axisym_lm_allocate_f_wav_multires(&f_wav_lm, &f_scal_lm, parameters);

  bandlimit = MIN(s2let_bandlimit(J_min - 1, parameters), L);
  ssht_adjoint_mw_inverse_sov_sym(
      f_scal_lm, f_scal, bandlimit, spin, dl_method, verbosity);
  offset = 0;
  offset_lm = 0;
  for (j = J_min; j <= J; j++) {
    int i;
    bandlimit = MIN(s2let_bandlimit(j, parameters), L);

    f_wav_norm = malloc(bandlimit * (2 * bandlimit - 1) * sizeof *f_wav_norm);

    for (i = 0; i < bandlimit * (2 * bandlimit - 1); ++i) {
      f_wav_norm[i] = f_wav[offset + i] / sqrt(2 * PI);
    }

    ssht_adjoint_mw_inverse_sov_sym(
        f_wav_lm + offset_lm, f_wav_norm, bandlimit, spin, dl_method, verbosity);
    offset_lm += bandlimit * bandlimit;
    offset += bandlimit * (2 * bandlimit - 1);

    free(f_wav_norm);
  }

  s2let_transform_axisym_lm_wav_synthesis_multires(
      flm, f_wav_lm, f_scal_lm, wav_lm, scal_lm, parameters);

  ssht_adjoint_mw_forward_sov_sym(f, flm, L, spin, dl_method, verbosity);

  free(flm);
  free(f_scal_lm);
  free(f_wav_lm);
}

/*!
 * Spherical wavelets : multiresolution synthesis in real space, MW sampling.
 * Perform multiresolution wavelet transform in real space (from scratch, gives wavelet
 * maps). Sampling scheme : MW sampling.
 *
 * \param[out]  f Input function (MW sampling)
 * \param[in]  f_wav Array of wavelets maps, MW sampling.
 * \param[in]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_wav_synthesis_mw_multires(
    complex double *f,
    const complex double *f_wav,
    const complex double *f_scal,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int spin = 0;
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;

  int bandlimit, j, offset, offset_lm;
  int J = s2let_j_max(parameters);
  // int l_min = s2let_transform_axisym_el_min(B, J_min);

  double *wav_lm, *scal_lm;
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, parameters);
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, parameters);

  complex double *flm, *f_wav_norm, *f_wav_lm, *f_scal_lm;
  flm = (complex double *)calloc(L * L, sizeof(complex double));
  s2let_transform_axisym_lm_allocate_f_wav_multires(&f_wav_lm, &f_scal_lm, parameters);

  bandlimit = MIN(s2let_bandlimit(J_min - 1, parameters), L);
  ssht_core_mw_forward_sov_conv_sym(
      f_scal_lm, f_scal, bandlimit, spin, dl_method, verbosity);
  offset = 0;
  offset_lm = 0;
  for (j = J_min; j <= J; j++) {
    int i;
    bandlimit = MIN(s2let_bandlimit(j, parameters), L);

    f_wav_norm = malloc(bandlimit * (2 * bandlimit - 1) * sizeof *f_wav_norm);

    for (i = 0; i < bandlimit * (2 * bandlimit - 1); ++i) {
      f_wav_norm[i] = f_wav[offset + i] * sqrt(2 * PI);
    }

    ssht_core_mw_forward_sov_conv_sym(
        f_wav_lm + offset_lm, f_wav_norm, bandlimit, spin, dl_method, verbosity);
    offset_lm += bandlimit * bandlimit;
    offset += bandlimit * (2 * bandlimit - 1);

    free(f_wav_norm);
  }

  s2let_transform_axisym_lm_wav_synthesis_multires(
      flm, f_wav_lm, f_scal_lm, wav_lm, scal_lm, parameters);

  ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

  free(flm);
  free(f_scal_lm);
  free(f_wav_lm);
}

/*!
 * Spherical wavelets : multiresolution synthesis adjoint in real space, MW sampling.
 * Perform multiresolution wavelet inverse transform adjoint in real space (from
 * scratch, gives wavelet maps). Sampling scheme : MW sampling.
 *
 * \param[out]  f_wav Array of wavelets maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  f Input function (MW sampling)
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_wav_synthesis_adjoint_mw_multires(
    complex double *f_wav,
    complex double *f_scal,
    const complex double *f,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int spin = 0;
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;

  int bandlimit, j, offset, offset_lm;
  int J = s2let_j_max(parameters);
  // int l_min = s2let_transform_axisym_el_min(B, J_min);

  double *wav_lm, *scal_lm;
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, parameters);
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, parameters);

  complex double *flm, *f_wav_lm, *f_scal_lm;
  flm = (complex double *)calloc(L * L, sizeof(complex double));
  s2let_transform_axisym_lm_allocate_f_wav_multires(&f_wav_lm, &f_scal_lm, parameters);

  ssht_adjoint_mw_inverse_sov_sym(flm, f, L, spin, dl_method, verbosity);

  s2let_transform_axisym_lm_wav_analysis_multires(
      f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, parameters);

  bandlimit = MIN(s2let_bandlimit(J_min - 1, parameters), L);
  ssht_adjoint_mw_forward_sov_sym(
      f_scal, f_scal_lm, bandlimit, spin, dl_method, verbosity);
  offset = 0;
  offset_lm = 0;
  for (j = J_min; j <= J; j++) {
    int i;

    bandlimit = MIN(s2let_bandlimit(j, parameters), L);
    ssht_adjoint_mw_forward_sov_sym(
        f_wav + offset, f_wav_lm + offset_lm, bandlimit, spin, dl_method, verbosity);

    for (i = 0; i < bandlimit * (2 * bandlimit - 1); ++i) {
      f_wav[offset + i] *= sqrt(2 * PI);
    }

    offset_lm += bandlimit * bandlimit;
    offset += bandlimit * (2 * bandlimit - 1);
  }

  free(flm);
  free(f_scal_lm);
  free(f_wav_lm);
}

/*!
 * Perform wavelet transform in real space (from scratch, gives wavelet maps).
 * Input function is real.
 * Sampling scheme : MW sampling.
 * Spherical wavelets : analysis in real space, MW sampling.
 *
 * \param[out]  f_wav Array of wavelets maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  f Input function (MW sampling)
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_wav_analysis_mw_real(
    double *f_wav,
    double *f_scal,
    const double *f,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;

  int j, offset, offset_lm;
  int J = s2let_j_max(parameters);
  // int l_min = s2let_transform_axisym_el_min(B, J_min);

  double *wav_lm, *scal_lm;
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, parameters);
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, parameters);

  complex double *flm, *f_wav_lm, *f_scal_lm;
  flm = (complex double *)calloc(L * L, sizeof(complex double));
  s2let_transform_axisym_lm_allocate_f_wav(&f_wav_lm, &f_scal_lm, parameters);

  ssht_core_mw_forward_sov_conv_sym_real(flm, f, L, dl_method, verbosity);

  s2let_transform_axisym_lm_wav_analysis(
      f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, parameters);

  ssht_core_mw_inverse_sov_sym_real(f_scal, f_scal_lm, L, dl_method, verbosity);
  offset = 0;
  offset_lm = 0;
  for (j = J_min; j <= J; j++) {
    int i;

    ssht_core_mw_inverse_sov_sym_real(
        f_wav + offset, f_wav_lm + offset_lm, L, dl_method, verbosity);

    for (i = 0; i < L * (2 * L - 1); ++i) {
      f_wav[offset + i] /= sqrt(2 * PI);
    }

    offset_lm += L * L;
    offset += L * (2 * L - 1);
  }

  free(flm);
  free(f_scal_lm);
  free(f_wav_lm);
}

/*!
 * Perform wavelet transform adjoint in real space (from scratch, gives wavelet maps).
 * Input function is real.
 * Sampling scheme : MW sampling.
 * Spherical wavelets : analysis adjoint in real space, MW sampling.
 *
 * \param[out]  f Input function (MW sampling)
 * \param[in]  f_wav Array of wavelets maps, MW sampling.
 * \param[in]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_wav_analysis_adjoint_mw_real(
    double *f,
    const double *f_wav,
    const double *f_scal,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;

  int j, offset, offset_lm;
  int J = s2let_j_max(parameters);
  // int l_min = s2let_transform_axisym_el_min(B, J_min);

  double *f_wav_norm, *wav_lm, *scal_lm;
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, parameters);
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, parameters);

  complex double *flm, *f_wav_lm, *f_scal_lm;
  flm = (complex double *)calloc(L * L, sizeof(complex double));
  s2let_transform_axisym_lm_allocate_f_wav(&f_wav_lm, &f_scal_lm, parameters);

  ssht_adjoint_mw_inverse_sov_sym_real(f_scal_lm, f_scal, L, dl_method, verbosity);
  offset = 0;
  offset_lm = 0;
  f_wav_norm = malloc(L * (2 * L - 1) * sizeof *f_wav_norm);
  for (j = J_min; j <= J; j++) {
    int i;

    for (i = 0; i < L * (2 * L - 1); ++i) {
      f_wav_norm[i] = f_wav[offset + i] / sqrt(2 * PI);
    }

    ssht_adjoint_mw_inverse_sov_sym_real(
        f_wav_lm + offset_lm, f_wav_norm, L, dl_method, verbosity);
    offset_lm += L * L;
    offset += L * (2 * L - 1);
  }
  free(f_wav_norm);

  s2let_transform_axisym_lm_wav_synthesis(
      flm, f_wav_lm, f_scal_lm, wav_lm, scal_lm, parameters);

  ssht_adjoint_mw_forward_sov_sym_real(f, flm, L, dl_method, verbosity);

  free(flm);
  free(f_scal_lm);
  free(f_wav_lm);
}

/*!
 * Perform wavelet transform in real space (from scratch, gives wavelet maps).
 * Input function is real.
 * Sampling scheme : MW sampling.
 * Spherical wavelets : synthesis in real space, MW sampling.
 *
 * \param[out]  f Input function (MW sampling)
 * \param[in]  f_wav Array of wavelets maps, MW sampling.
 * \param[in]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_wav_synthesis_mw_real(
    double *f,
    const double *f_wav,
    const double *f_scal,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;

  int j, offset, offset_lm;
  int J = s2let_j_max(parameters);
  // int l_min = s2let_transform_axisym_el_min(B, J_min);

  double *f_wav_norm, *wav_lm, *scal_lm;
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, parameters);
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, parameters);

  complex double *flm, *f_wav_lm, *f_scal_lm;
  flm = (complex double *)calloc(L * L, sizeof(complex double));
  s2let_transform_axisym_lm_allocate_f_wav(&f_wav_lm, &f_scal_lm, parameters);

  ssht_core_mw_forward_sov_conv_sym_real(f_scal_lm, f_scal, L, dl_method, verbosity);
  offset = 0;
  offset_lm = 0;
  f_wav_norm = malloc(L * (2 * L - 1) * sizeof *f_wav_norm);
  for (j = J_min; j <= J; j++) {
    int i;

    for (i = 0; i < L * (2 * L - 1); ++i) {
      f_wav_norm[i] = f_wav[offset + i] * sqrt(2 * PI);
    }

    ssht_core_mw_forward_sov_conv_sym_real(
        f_wav_lm + offset_lm, f_wav_norm, L, dl_method, verbosity);
    offset_lm += L * L;
    offset += L * (2 * L - 1);
  }
  free(f_wav_norm);

  s2let_transform_axisym_lm_wav_synthesis(
      flm, f_wav_lm, f_scal_lm, wav_lm, scal_lm, parameters);

  ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);

  free(flm);
  free(f_scal_lm);
  free(f_wav_lm);
}

/*!
 * Perform wavelet inverse transform adjoint in real space (from scratch, gives wavelet
 * maps). Input function is real. Sampling scheme : MW sampling. Spherical wavelets :
 * analysis in real space, MW sampling.
 *
 * \param[out]  f_wav Array of wavelets maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  f Input function (MW sampling)
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_wav_synthesis_adjoint_mw_real(
    double *f_wav,
    double *f_scal,
    const double *f,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;

  int j, offset, offset_lm;
  int J = s2let_j_max(parameters);
  // int l_min = s2let_transform_axisym_el_min(B, J_min);

  double *wav_lm, *scal_lm;
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, parameters);
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, parameters);

  complex double *flm, *f_wav_lm, *f_scal_lm;
  flm = (complex double *)calloc(L * L, sizeof(complex double));
  s2let_transform_axisym_lm_allocate_f_wav(&f_wav_lm, &f_scal_lm, parameters);

  ssht_adjoint_mw_inverse_sov_sym_real(flm, f, L, dl_method, verbosity);

  s2let_transform_axisym_lm_wav_analysis(
      f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, parameters);

  ssht_adjoint_mw_forward_sov_sym_real(f_scal, f_scal_lm, L, dl_method, verbosity);
  offset = 0;
  offset_lm = 0;
  for (j = J_min; j <= J; j++) {
    int i;

    ssht_adjoint_mw_forward_sov_sym_real(
        f_wav + offset, f_wav_lm + offset_lm, L, dl_method, verbosity);

    for (i = 0; i < L * (2 * L - 1); ++i) {
      f_wav[offset + i] *= sqrt(2 * PI);
    }

    offset_lm += L * L;
    offset += L * (2 * L - 1);
  }

  free(flm);
  free(f_scal_lm);
  free(f_wav_lm);
}

/*!
 * Spherical wavelets : multiresolution analysis in real space, MW sampling.
 * Perform multiresolution wavelet transform in real space (from scratch, gives wavelet
 * maps). Input function is real. Sampling scheme : MW sampling.
 *
 * \param[out]  f_wav Array of wavelets maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  f Input function (MW sampling)
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_wav_analysis_mw_multires_real(
    double *f_wav,
    double *f_scal,
    const double *f,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;

  int bandlimit, j, offset, offset_lm;
  int J = s2let_j_max(parameters);
  // int l_min = s2let_transform_axisym_el_min(B, J_min);

  double *wav_lm, *scal_lm;
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, parameters);
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, parameters);

  complex double *flm, *f_wav_lm, *f_scal_lm;
  flm = (complex double *)calloc(L * L, sizeof(complex double));
  s2let_transform_axisym_lm_allocate_f_wav_multires(&f_wav_lm, &f_scal_lm, parameters);

  ssht_core_mw_forward_sov_conv_sym_real(flm, f, L, dl_method, verbosity);

  s2let_transform_axisym_lm_wav_analysis_multires(
      f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, parameters);

  bandlimit = MIN(s2let_bandlimit(J_min - 1, parameters), L);
  ssht_core_mw_inverse_sov_sym_real(f_scal, f_scal_lm, bandlimit, dl_method, verbosity);
  offset = 0;
  offset_lm = 0;
  for (j = J_min; j <= J; j++) {
    int i;

    bandlimit = MIN(s2let_bandlimit(j, parameters), L);
    ssht_core_mw_inverse_sov_sym_real(
        f_wav + offset, f_wav_lm + offset_lm, bandlimit, dl_method, verbosity);

    for (i = 0; i < bandlimit * (2 * bandlimit - 1); ++i) {
      f_wav[offset + i] /= sqrt(2 * PI);
    }

    offset_lm += bandlimit * bandlimit;
    offset += bandlimit * (2 * bandlimit - 1);
  }

  free(flm);
  free(f_scal_lm);
  free(f_wav_lm);
}

/*!
 * Spherical wavelets : multiresolution analysis adjoint in real space, MW sampling.
 * Perform wavelet transform adjoint in real space (from scratch, gives wavelet maps).
 * Input function is real.
 * Sampling scheme : MW sampling.
 *
 * \param[out]  f Input function (MW sampling)
 * \param[in]  f_wav Array of wavelets maps, MW sampling.
 * \param[in]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_wav_analysis_adjoint_mw_multires_real(
    double *f,
    const double *f_wav,
    const double *f_scal,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;

  int bandlimit, j, offset, offset_lm;
  int J = s2let_j_max(parameters);
  // int l_min = s2let_transform_axisym_el_min(B, J_min);

  double *f_wav_norm, *wav_lm, *scal_lm;
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, parameters);
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, parameters);

  complex double *flm, *f_wav_lm, *f_scal_lm;
  flm = (complex double *)calloc(L * L, sizeof(complex double));
  s2let_transform_axisym_lm_allocate_f_wav_multires(&f_wav_lm, &f_scal_lm, parameters);

  bandlimit = MIN(s2let_bandlimit(J_min - 1, parameters), L);
  ssht_adjoint_mw_inverse_sov_sym_real(
      f_scal_lm, f_scal, bandlimit, dl_method, verbosity);
  offset = 0;
  offset_lm = 0;
  for (j = J_min; j <= J; j++) {
    int i;
    bandlimit = MIN(s2let_bandlimit(j, parameters), L);

    f_wav_norm = malloc(bandlimit * (2 * bandlimit - 1) * sizeof *f_wav_norm);

    for (i = 0; i < bandlimit * (2 * bandlimit - 1); ++i) {
      f_wav_norm[i] = f_wav[offset + i] / sqrt(2 * PI);
    }

    ssht_adjoint_mw_inverse_sov_sym_real(
        f_wav_lm + offset_lm, f_wav_norm, bandlimit, dl_method, verbosity);
    offset_lm += bandlimit * bandlimit;
    offset += bandlimit * (2 * bandlimit - 1);

    free(f_wav_norm);
  }

  s2let_transform_axisym_lm_wav_synthesis_multires(
      flm, f_wav_lm, f_scal_lm, wav_lm, scal_lm, parameters);

  ssht_adjoint_mw_forward_sov_sym_real(f, flm, L, dl_method, verbosity);

  free(flm);
  free(f_scal_lm);
  free(f_wav_lm);
}

/*!
 * Spherical wavelets : multiresolution synthesis in real space, MW sampling.
 * Perform wavelet transform in real space (from scratch, gives wavelet maps).
 * Input function is real.
 * Sampling scheme : MW sampling.
 *
 * \param[out]  f Input function (MW sampling)
 * \param[in]  f_wav Array of wavelets maps, MW sampling.
 * \param[in]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_wav_synthesis_mw_multires_real(
    double *f,
    const double *f_wav,
    const double *f_scal,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;

  int bandlimit, j, offset, offset_lm;
  int J = s2let_j_max(parameters);
  // int l_min = s2let_transform_axisym_el_min(B, J_min);

  double *f_wav_norm, *wav_lm, *scal_lm;
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, parameters);
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, parameters);

  complex double *flm, *f_wav_lm, *f_scal_lm;
  flm = (complex double *)calloc(L * L, sizeof(complex double));
  s2let_transform_axisym_lm_allocate_f_wav_multires(&f_wav_lm, &f_scal_lm, parameters);

  bandlimit = MIN(s2let_bandlimit(J_min - 1, parameters), L);
  ssht_core_mw_forward_sov_conv_sym_real(
      f_scal_lm, f_scal, bandlimit, dl_method, verbosity);
  offset = 0;
  offset_lm = 0;
  for (j = J_min; j <= J; j++) {
    int i;
    bandlimit = MIN(s2let_bandlimit(j, parameters), L);

    f_wav_norm = malloc(bandlimit * (2 * bandlimit - 1) * sizeof *f_wav_norm);

    for (i = 0; i < bandlimit * (2 * bandlimit - 1); ++i) {
      f_wav_norm[i] = f_wav[offset + i] * sqrt(2 * PI);
    }

    ssht_core_mw_forward_sov_conv_sym_real(
        f_wav_lm + offset_lm, f_wav_norm, bandlimit, dl_method, verbosity);
    offset_lm += bandlimit * bandlimit;
    offset += bandlimit * (2 * bandlimit - 1);

    free(f_wav_norm);
  }

  s2let_transform_axisym_lm_wav_synthesis_multires(
      flm, f_wav_lm, f_scal_lm, wav_lm, scal_lm, parameters);

  ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);

  free(flm);
  free(f_scal_lm);
  free(f_wav_lm);
}

/*!
 * Spherical wavelets : multiresolution synthesis adjoint in real space, MW sampling.
 * Perform multiresolution wavelet inverse transform adjoint in real space (from
 * scratch, gives wavelet maps). Input function is real. Sampling scheme : MW sampling.
 *
 * \param[out]  f_wav Array of wavelets maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  f Input function (MW sampling)
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_wav_synthesis_adjoint_mw_multires_real(
    double *f_wav,
    double *f_scal,
    const double *f,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;

  int bandlimit, j, offset, offset_lm;
  int J = s2let_j_max(parameters);
  // int l_min = s2let_transform_axisym_el_min(B, J_min);

  double *wav_lm, *scal_lm;
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, parameters);
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, parameters);

  complex double *flm, *f_wav_lm, *f_scal_lm;
  flm = (complex double *)calloc(L * L, sizeof(complex double));
  s2let_transform_axisym_lm_allocate_f_wav_multires(&f_wav_lm, &f_scal_lm, parameters);

  ssht_adjoint_mw_inverse_sov_sym_real(flm, f, L, dl_method, verbosity);

  s2let_transform_axisym_lm_wav_analysis_multires(
      f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, parameters);

  bandlimit = MIN(s2let_bandlimit(J_min - 1, parameters), L);
  ssht_adjoint_mw_forward_sov_sym_real(
      f_scal, f_scal_lm, bandlimit, dl_method, verbosity);
  offset = 0;
  offset_lm = 0;
  for (j = J_min; j <= J; j++) {
    int i;

    bandlimit = MIN(s2let_bandlimit(j, parameters), L);
    ssht_adjoint_mw_forward_sov_sym_real(
        f_wav + offset, f_wav_lm + offset_lm, bandlimit, dl_method, verbosity);

    for (i = 0; i < bandlimit * (2 * bandlimit - 1); ++i) {
      f_wav[offset + i] *= sqrt(2 * PI);
    }

    offset_lm += bandlimit * bandlimit;
    offset += bandlimit * (2 * bandlimit - 1);
  }

  free(flm);
  free(f_scal_lm);
  free(f_wav_lm);
}

/*!
 * Threshold real wavelets in real space, MW sampling, full resolution.
 *
 * \param[inout]  g_wav Array of wavelets maps, MW sampling.
 * \param[in]  threshold A threshold rule, i.e. a number for every scale j.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_wav_hardthreshold_real(
    double *g_wav, const double *threshold, const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int J = s2let_j_max(parameters);
  int i, j, offset = 0;
  for (j = J_min; j <= J; j++) {
    int bl = L;
    for (i = 0; i < bl * (2 * bl - 1); i++) {
      if (abs(g_wav[offset + i]) < threshold[j - J_min])
        g_wav[offset + i] = 0;
    }
    offset += bl * (2 * bl - 1);
  }
}

/*!
 * Threshold real wavelets in real space, MW sampling, multiresolution.
 *
 * \param[inout]  g_wav Array of wavelets maps, MW sampling.
 * \param[in]  threshold A threshold rule, i.e. a number for every scale j.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_wav_hardthreshold_multires_real(
    double *g_wav, const double *threshold, const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int J = s2let_j_max(parameters);
  int i, j, offset = 0;
  for (j = J_min; j <= J; j++) {
    int bl = MIN(s2let_bandlimit(j, parameters), L);
    for (i = 0; i < bl * (2 * bl - 1); i++) {
      if (abs(g_wav[offset + i]) < threshold[j - J_min])
        g_wav[offset + i] = 0;
    }
    offset += bl * (2 * bl - 1);
  }
}
