// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_AXISYM_DENOISE_DEMO
#define S2LET_AXISYM_DENOISE_DEMO
#include <ssht/ssht.h>

#ifdef __cplusplus
extern "C" {
#endif

//!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void s2let_lm_random_flm_real_sigma(complex double *flm, int L, int seed, double sigmanoise);

//!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double needletpower(double *wav_lm, int L);


#ifdef __cplusplus
}
#endif
#endif
