// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_DENOISE_DEMO
#define S2LET_DENOISE_DEMO
#include <ssht/ssht.h>

#ifdef __cplusplus
extern "C" {
#endif

//!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void s2let_lm_random_flm_real_sigma(complex double *flm, int L, int seed, double sigmanoise);

//!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double waveletpower(complex double *wav_lm, int L);

//!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void hard_threshold_real(double *g_wav, const double *threshold, const s2let_parameters_t *parameters); 

#ifdef __cplusplus
}
#endif
#endif
