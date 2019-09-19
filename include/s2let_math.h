// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_MATH
#define S2LET_MATH
#include <ssht/ssht.h>

#define PI 3.141592653589793238462643383279502884197

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define ABS(a) ((a) > (0.0) ? (a) : (-(a)))

#ifdef __cplusplus
extern "C" {
#endif

double s2let_math_kappa0_quadtrap_s2dw(double a, double b, int n, double B);

double s2let_math_kappa0_quadtrap_needlet(double a, double b, int n);

double s2let_math_spline_scalingfct(double x, double y);

double ran2_dp(int idum);

double maxerr_cplx(S2LET_COMPLEX(double) * a, S2LET_COMPLEX(double) * b, int size);

double maxerr(double* a, double* b, int size);

unsigned long binomial_coefficient(int n, int k, int exact);

#ifdef __cplusplus
}
#endif

#endif
