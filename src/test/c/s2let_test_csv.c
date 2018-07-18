// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Büttner and Jason McEwen
// See LICENSE.txt for license details

/*!
 * \file so3_test_csv.c
 * Applies SO3 algorithms to perform inverse and forward Wigner
 * transforms (respectively) to check that the original
 * signal is reconstructed exactly (to numerical precision). Test is
 * performed on a random signal with harmonic coefficients uniformly
 * sampled from (-1,1), once for a real and once for a complex signal.
 * This test is run for multiple values of L and outputs the results
 * in CSV format (to stdout).
 * If no value for N is given in the arguments, each test will run
 * with N = L, otherwise all tests will be run with the given N (and
 * tests with L < N will be skipped). L will take all powers of 2 less
 * or equal than the given Lmax. If L0 is given, all powers of 2 less
 * than or equal to L0 will be skipped (and L0 will be used as the lower
 * band-limit for all transforms).
 *
 * \par Usage
 *   \code{.sh}
 *   so3_test_csv [Lmax [L0 [N]]]
 *   \endcode
 *   e.g.
 *   \code{.sh}
 *   so3_test 128 0 4
 *   \endcode
 *   Defaults: L = 64, L0 = 0, N = L<br>
 *   To save it to a file, just use something like
 *   \code{.sh}
 *   so3_test_csv 64 7 4 > test.csv
 *   \endcode
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>

#include <ssht.h>
#include <s2let.h>

#define NREPEAT 10

double get_max_error(complex double *expected, complex double *actual, int n);
double ran2_dp(int idum);

int main(int argc, char **argv)
{
    s2let_parameters_t parameters = {};
    int verbosity = parameters.verbosity = 0;
    int Lmax, L, useLasN, N, L0, B, J_min, spin;
    complex *flm, *flm_rec, *f, *f_rec;;
    complex *f_wav, *f_scal;
    int seed;
    clock_t time_start, time_end;
    int i, multires;

    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;

    double min_duration_inverse;
    double min_duration_forward;
    double avg_duration_inverse;
    double avg_duration_forward;
    double avg_error;

    // Parse command line arguments
    parameters.N = N = Lmax = 64;
    useLasN = 1; // true
    L0 = 0;
    parameters.spin = spin = 0;
    parameters.J_min = J_min = 0;
    seed = 1;
    parameters.B = B = 2;
    if (argc > 1)
        Lmax = atoi(argv[1]);
    if (argc > 2)
        parameters.spin = spin = atoi(argv[2]);
    if (argc > 3)
        parameters.B = B = atoi(argv[3]);
    if (argc > 4)
    {
        useLasN = 0; // false
        parameters.N = N = atoi(argv[4]);
    }

    parameters.original_spin = 0;

    // Output header row
    printf("multires;spin;L;N;B;J_min;J;min_duration_inverse;min_duration_forward;avg_duration_inverse;avg_duration_forward;avg_error\n");

    // multires == 0 --> full resolution transform
    // multires == 1 --> multiresolution transform
    for (multires = 0; multires < 2; ++multires)
    {
        parameters.upsample = 1-multires;
        L = 1024;
        while(L <= Lmax)
        {
            if (L <= L0 || (!useLasN && L < N))
            {
                L *= 2;
                continue;
            }

            if (useLasN)
            {
                parameters.N = N = L;
            }

            parameters.L = L;

            s2let_allocate_lm(&flm, L);
            s2let_allocate_lm(&flm_rec, L);
            s2let_allocate_mw(&f, L);
            s2let_allocate_mw(&f_rec, L);
            s2let_allocate_f_wav(&f_wav, &f_scal, &parameters);

            min_duration_inverse = 0.0;
            min_duration_forward = 0.0;
            avg_duration_inverse = 0.0;
            avg_duration_forward = 0.0;
            avg_error = 0.0;
            int J = s2let_j_max(&parameters);

            for (i = 0; i < NREPEAT; ++i)
            {
                int j;
                double duration;

                // Reset output array
                for (j = 0; j < L*L; ++j)
                    flm_rec[j] = 0.0;

        s2let_lm_random_flm(flm, L, spin, seed);
        ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

                time_start = clock();
                s2let_analysis_px2wav(f_wav, f_scal, f, &parameters);
                time_end = clock();

                duration = (time_end - time_start) / (double)CLOCKS_PER_SEC;
		avg_duration_inverse += duration / NREPEAT;
                if (!i || duration < min_duration_inverse)
                    min_duration_inverse = duration;

                time_start = clock();
                s2let_synthesis_wav2px(f_rec, f_wav, f_scal, &parameters);
                time_end = clock();

        ssht_core_mw_forward_sov_conv_sym(flm_rec, f_rec, L, spin, dl_method, verbosity);

                duration = (time_end - time_start) / (double)CLOCKS_PER_SEC;
		avg_duration_forward += duration / NREPEAT;
                if (!i || duration < min_duration_forward)
                    min_duration_forward = duration;

                avg_error += get_max_error(flm, flm_rec, L*L)/NREPEAT;
            }

            printf("%d;%d;%d;%d;%d;%d;%d;%f;%f;%f;%f;%e\n",
                   multires,
           spin,
                   L,
           N,
           B,
                   J_min,
                   J,
                   min_duration_inverse,
                   min_duration_forward,
                   avg_duration_inverse,
                   avg_duration_forward,
                   avg_error);

            L *= 2;

        free(f);
        free(f_rec);
        free(flm);
        free(flm_rec);
        free(f_wav);
        free(f_scal);
        }

    }


    return 0;
}

double get_max_error(complex double *expected, complex double *actual, int n)
{
    int i;
    double error, maxError = 0;

    for (i = 0; i < n; ++i)
    {
        error = cabs(expected[i] - actual[i]);
        maxError = MAX(error, maxError);
    }

    return maxError;
}
