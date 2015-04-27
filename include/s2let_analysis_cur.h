// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_ANALYSIS_CUR   //
#define S2LET_ANALYSIS_CUR   //

#include <complex.h>

/** Harmonic-space curvelet transform **/

void s2let_analysis_cur_lm2lmn(
    complex double *f_cur_lmn,
    complex double *f_scal_lm,
    const complex double *flm,
    const complex double *cur_lm,
    const double *scal_l,
    const s2let_parameters_t *parameters
);

/** Harmonic-space curvelet transform for real signals **/

void s2let_analysis_cur_lm2lmn_real(
    complex double *f_cur_lmn,
    complex double *f_scal_lm,
    const complex double *flm,
    const complex double *cur_lm,
    const double *scal_l,
    const s2let_parameters_t *parameters
);

/** Harmonic-space curvelet transform **/

void s2let_analysis_lm2cur(
    complex double *f_cur,
    complex double *f_scal,
    const complex double *flm,
    const s2let_parameters_t *parameters
);

/** Harmonic-space curvelet transform for real signals **/

void s2let_analysis_lm2cur_real(
    double *f_cur,
    double *f_scal,
    const complex double *flm,
    const s2let_parameters_t *parameters
);

/** Pixel-space curvelet transform **/

void s2let_analysis_px2cur(
    complex double *f_cur,
    complex double *f_scal,
    const complex double *f,
    const s2let_parameters_t *parameters
);

/** Pixel-space curvelet transform for real signals**/

void s2let_analysis_px2cur_real(
    double *f_cur,
    double *f_scal,
    const double *f,
    const s2let_parameters_t *parameters
);

#endif
