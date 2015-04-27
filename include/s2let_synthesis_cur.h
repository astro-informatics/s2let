// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_SYNTHESIS_CUR      //
#define S2LET_SYNTHESIS_CUR       // 

#include <complex.h>

/** Harmonic-space curvelet transform **/

void s2let_synthesis_cur_lmn2lm(
    complex double *flm,
    const complex double *f_cur_lmn,
    const complex double *f_scal_lm,
    const complex double *cur_lm,
    const double *scal_l,
    const s2let_parameters_t *parameters
);

/** Harmonic-space curvelet transform for real signals **/

void s2let_synthesis_cur_lmn2lm_real(
    complex double *flm,
    const complex double *f_cur_lmn,
    const complex double *f_scal_lm,
    const complex double *cur_lm,
    const double *scal_l,
    const s2let_parameters_t *parameters
);

/** Harmonic-space curvelet transform **/

void s2let_synthesis_cur2lm(
    complex double *flm,
    const complex double *f_cur,
    const complex double *f_scal,
    const s2let_parameters_t *parameters
);

/** Harmonic-space curvelet transform for real signals **/

void s2let_synthesis_cur2lm_real(
    complex double *flm,
    const double *f_cur,
    const double *f_scal,
    const s2let_parameters_t *parameters
);

/** Pixel-space curvelet transform **/

void s2let_synthesis_cur2px(
    complex double *f,
    const complex double *f_cur,
    const complex double *f_scal,
    const s2let_parameters_t *parameters
);

/** Pixel-space curvelet transform for real signals**/

void s2let_synthesis_cur2px_real(
    double *f,
    const double *f_cur,
    const double *f_scal,
    const s2let_parameters_t *parameters
);

#endif
