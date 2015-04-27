// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_ALLOC
#define S2LET_ALLOC

#include <complex.h>

/** Pixel space allocation **/

void s2let_allocate_mw(complex double **f, int L);
void s2let_allocate_mw_real(double **f, int L);

void s2let_allocate_mwss(complex double **f, int L);
void s2let_allocate_mwss_real(double **f, int L);

/** Harmonic space allocation **/

void s2let_allocate_lm(complex double **flm, int L);

/** Wigner space allocation **/

void s2let_allocate_lmn_f_wav(
    complex double **f_wav_lmn,
    complex double **f_scal_lm,
    const s2let_parameters_t *parameters
);

/** Wavelet space allocation **/

void s2let_allocate_f_wav(
    complex double **f_wav,
    complex double **f_scal,
    const s2let_parameters_t *parameters
);

void s2let_allocate_f_wav_real(
    double **f_wav,
    double **f_scal,
    const s2let_parameters_t *parameters
);



/** Wigner space allocation **/
void s2let_allocate_lmn_f_cur(
    complex double **f_cur_lmn,
    complex double **f_scal_lm,
    const s2let_parameters_t *parameters
);

/** Curvelet space allocation **/
void s2let_allocate_f_cur(
     complex double **f_cur,
     complex double **f_scal,
     const s2let_parameters_t *parameters
);

void s2let_allocate_f_cur_real(
    double **f_cur,
    double **f_scal,
    const s2let_parameters_t *parameters
);



#endif
