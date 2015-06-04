// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_MW
#define S2LET_MW


/** Interfaces to SSHT (required by the Java interface to S2LET) **/
void s2let_mw_map2alm_real(complex double* flm, const double* f, int L);
void s2let_mw_alm2map_real(double* f, const complex double* flm, int L);
void s2let_mw_alm2map(complex double* f, const complex double* flm, int L, int spin);
void s2let_mw_map2alm(complex double* flm, const complex double* f, int L, int spin);

/** Helper functions for pixel-space computations in MW sampling **/
double s2let_mw_power(complex double *flm, int L);
double s2let_mw_power_real(double *flm, int L);

#endif
