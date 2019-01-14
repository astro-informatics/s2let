// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <complex.h>
#include <ssht.h>
#include <so3.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*!
 * Compute power for a signal in harmonic space.
 */
double s2let_lm_power(complex double *flm, int L){
    int i;
    double totalpower = 0.0;
    for(i = 0; i < L*L; ++i)
        totalpower += pow(cabs(flm[i]), 2.0);
    totalpower = totalpower / (L * L);
    return totalpower;
}

/*!
 * Generate random harmonic coefficients for a complex map.
 *
 * \param[out]  flm Harmonic coefficients.
 * \param[in]  L Band-limit.
 * \parma[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_lm_random_flm(complex double *flm, int L, int spin, int seed)
{
    int i, i_min;
    srand( time(NULL) );
    // el < |s| are zero, so start with el = |s|, m = -el.
    i_min = spin*spin;
    for (i=i_min; i<L*L; ++i)
        flm[i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
}

/*!
 * Generate random harmonic coefficients corresponding to a real map.
 *
 * \param[out]  flm Harmonic coefficients.
 * \param[in]  L Band-limit.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_lm_random_flm_real(complex double *flm, int L, int seed) {
    int el, m, msign, i, i_op;
    for (el = 0; el < L; ++el) {
        m = 0;
        i = el*el + el + m;
        flm[i] = (2.0*ran2_dp(seed) - 1.0);
        for (m = 1; m <= el; ++m) {
            i = el*el + el + m;
            flm[i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
            i_op = el*el + el - m;
            msign = m & 1;
            msign = 1 - msign - msign; // (-1)^m
            flm[i_op] = msign * conj(flm[i]);
        }
    }
}
