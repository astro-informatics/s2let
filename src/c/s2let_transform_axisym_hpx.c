// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include <complex.h>
#include <stdlib.h>

#include "s2let/s2let.h"

void s2let_transform_axisym_allocate_hpx_f_wav_real(
    double **f_wav, double **f_scal, int nside, const s2let_parameters_t *parameters) {
  int J_min = parameters->J_min;

  int J = s2let_j_max(parameters);
  *f_wav = (double *)calloc((J + 1 - J_min) * 12 * nside * nside, sizeof(double));
  *f_scal = (double *)calloc(12 * nside * nside, sizeof(double));
}

void s2let_transform_axisym_wav_analysis_hpx_real(
    double *f_wav,
    double *f_scal,
    const double *f,
    int nside,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int bandlimit, j, offset, offset_lm;
  int J = s2let_j_max(parameters);

  double *wav_lm, *scal_lm;
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, parameters);
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, parameters);

  complex double *flm, *f_wav_lm, *f_scal_lm;
  flm = (complex double *)calloc(L * L, sizeof(complex double));
  s2let_transform_axisym_lm_allocate_f_wav_multires(&f_wav_lm, &f_scal_lm, parameters);

  s2let_hpx_map2alm_real(flm, f, nside, L);

  s2let_transform_axisym_lm_wav_analysis_multires(
      f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, parameters);

  bandlimit = MIN(s2let_bandlimit(J_min - 1, parameters), L);
  s2let_hpx_alm2map_real(f_scal, f_scal_lm, nside, bandlimit);

  offset = 0;
  offset_lm = 0;
  for (j = J_min; j <= J; j++) {
    bandlimit = MIN(s2let_bandlimit(j, parameters), L);
    s2let_hpx_alm2map_real(f_wav + offset, f_wav_lm + offset_lm, nside, bandlimit);
    offset_lm += bandlimit * bandlimit;
    offset += 12 * nside * nside;
  }

  free(flm);
  free(f_scal_lm);
  free(f_wav_lm);
}

void s2let_transform_axisym_wav_synthesis_hpx_real(
    double *f,
    const double *f_wav,
    const double *f_scal,
    int nside,
    const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;

  int bandlimit, j, offset, offset_lm;
  int J = s2let_j_max(parameters);

  double *wav_lm, *scal_lm;
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, parameters);
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, parameters);

  complex double *flm, *f_wav_lm, *f_scal_lm;
  flm = (complex double *)calloc(L * L, sizeof(complex double));
  s2let_transform_axisym_lm_allocate_f_wav_multires(&f_wav_lm, &f_scal_lm, parameters);

  bandlimit = MIN(s2let_bandlimit(J_min - 1, parameters), L);
  s2let_hpx_map2alm_real(f_scal_lm, f_scal, nside, bandlimit);
  offset = 0;
  offset_lm = 0;
  for (j = J_min; j <= J; j++) {
    bandlimit = MIN(s2let_bandlimit(j, parameters), L);
    s2let_hpx_map2alm_real(f_wav_lm + offset_lm, f_wav + offset, nside, bandlimit);
    offset_lm += bandlimit * bandlimit;
    offset += 12 * nside * nside;
  }

  s2let_transform_axisym_lm_wav_synthesis_multires(
      flm, f_wav_lm, f_scal_lm, wav_lm, scal_lm, parameters);

  s2let_hpx_alm2map_real(f, flm, nside, L);

  free(flm);
  free(f_scal_lm);
  free(f_wav_lm);
}
