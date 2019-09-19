// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_MW
#define S2LET_MW

#include <ssht/ssht.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Interfaces to SSHT (required by the Java interface to S2LET) **/
void s2let_mw_map2alm_real(S2LET_COMPLEX(double) * flm, const double* f, int L);
void s2let_mw_alm2map_real(double* f, const S2LET_COMPLEX(double) * flm, int L);
void s2let_mw_alm2map(S2LET_COMPLEX(double) * f, const S2LET_COMPLEX(double) * flm, int L, int spin);
void s2let_mw_map2alm(S2LET_COMPLEX(double) * flm, const S2LET_COMPLEX(double) * f, int L, int spin);

/** Helper functions for pixel-space computations in MW sampling **/
double s2let_mw_power(S2LET_COMPLEX(double) * flm, int L);
double s2let_mw_power_real(double* flm, int L);

#ifdef __cplusplus
}
#endif
#endif
