// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_LM
#define S2LET_LM

#include <ssht/ssht.h>

#ifdef __cplusplus
extern "C" {
#endif

double s2let_lm_power(S2LET_COMPLEX(double) *flm, int L);

void s2let_lm_random_flm(S2LET_COMPLEX(double) *flm, int L, int spin, int seed);
void s2let_lm_random_flm_real(S2LET_COMPLEX(double) *flm, int L, int seed);

#ifdef __cplusplus
}
#endif
#endif
