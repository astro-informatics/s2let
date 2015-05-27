// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <ssht.h>
#include <assert.h>
#include <complex.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*!
 * Test the identity relation of the wavelet tiling in harmonic space.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_tiling_axisym_test(int B, int L, int J_min)
{
  s2let_parameters_t parameters = {};
  parameters.B = B;
  parameters.L = L;
  parameters.J_min = J_min;

  double *kappa, *kappa0;

  // Allocate the kernels corresponding to the parameters B, L
  s2let_tiling_axisym_allocate(&kappa, &kappa0, &parameters);

  // Construct the tiling of harmonic space
  s2let_tiling_axisym(kappa, kappa0, &parameters);

  // Check that they recover the identity relation,
  // ensuring exactness of the wavelet transform.
  double res = s2let_tiling_axisym_check_identity(kappa, kappa0, &parameters);
  printf("  - Maximum error : %6.5e\n", res);

  free(kappa);
  free(kappa0);
}



/*!
 * Test the identity relation of the directionality components for
 * the wavelet tiling in harmonic space.
 *
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  N Azimuthal band-limit.
 * \retval none
 */
void s2let_tiling_curvelet_direction_test(int L, int N)
{
    s2let_parameters_t parameters = {};
    parameters.L = L;
    parameters.N = N;
    
    complex double *s_elm;
    double error;
    
    // Allocate space for the harmonic coefficients
    s2let_tiling_direction_allocate(&s_elm, &parameters);
    
    // Construct the harmonic coefficients
    s2let_tiling_curvelet_direction(s_elm, &parameters);
    
    // Check that they recover the identity relation,
    // ensuring exactness of the wavelet transform.
    error = s2let_tiling_direction_check_identity(s_elm, &parameters);
    printf("  - Maximum error : %6.5e\n", error);
    
    free(s_elm);
}


/*!
 * Test the identity relation of the directionality components for
 * the wavelet tiling in harmonic space.
 *
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  N Azimuthal band-limit.
 * \retval none
 */
void s2let_tiling_direction_test(int L, int N)
{
  s2let_parameters_t parameters = {};
  parameters.L = L;
  parameters.N = N;

  complex double *s_elm;
  double error;

  // Allocate space for the harmonic coefficients
  s2let_tiling_direction_allocate(&s_elm, &parameters);

  // Construct the harmonic coefficients
  s2let_tiling_direction(s_elm, &parameters);

  // Check that they recover the identity relation,
  // ensuring exactness of the wavelet transform.
  error = s2let_tiling_direction_check_identity(s_elm, &parameters);
  printf("  - Maximum error : %6.5e\n", error);

  free(s_elm);
}

/*!
 * Test the identity relation of the directional wavelets for
 * the wavelet tiling in harmonic space.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \retval none
 */
void s2let_tiling_wavelet_test(int B, int L, int J_min, int N, int spin)
{
  s2let_parameters_t parameters = {};
  parameters.B = B;
  parameters.L = L;
  parameters.J_min = J_min;
  parameters.N = N;
  parameters.spin = spin;
  parameters.normalization = S2LET_WAV_NORM_DEFAULT;
  parameters.original_spin = 0;

  complex double *phi;
  double *psi;
  double error;

  // Allocate space for the harmonic coefficients
  s2let_tiling_wavelet_allocate(&phi, &psi, &parameters);

  // Construct the harmonic coefficients
  s2let_tiling_wavelet(phi, psi, &parameters);

  // Check that they recover the identity relation,
  // ensuring exactness of the wavelet transform.
  error = s2let_tiling_wavelet_check_identity(phi, psi, &parameters);
  printf("  - Maximum error : %6.5e\n", error);

  free(phi);
  free(psi);
}



/*!
 * Test the identity relation of the curvelets for
 * the wavelet tiling in harmonic space.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \retval none
 */
void s2let_tiling_curvelet_test(int B, int L, int J_min, int N, int spin)
{
    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = spin;
    parameters.normalization = S2LET_WAV_NORM_DEFAULT;
    parameters.original_spin = 0;
    
    complex double *phi;
    double *psi;
    double error;
    
    // Allocate space for the harmonic coefficients
    s2let_tiling_curvelet_allocate(&phi, &psi, &parameters);
    
    // Construct the harmonic coefficients
    s2let_tiling_curvelet(phi, psi, &parameters);
    
    // Check that they recover the identity relation,
    // ensuring exactness of the wavelet transform.
    error = s2let_tiling_curvelet_check_identity(phi, psi, &parameters);
    printf("  - Maximum error : %6.5e\n", error);
    
    free(phi);
    free(psi);
}



/*!
 * Test binomial coefficient
 */

void s2let_binomial_coefficient_test(int n_max)
{
    const int nRepeat = 100000;
    int n, k, i;
    int firstError = 0;

    long error;

    clock_t time_start, time_end;

    for (n = 1; n <= n_max; ++n)
    {
        for (k = 0; k <= n/2; ++k)
        {
            error = binomial_coefficient(n, k, 0) -
                        binomial_coefficient(n, k, 1);

            if (error && !firstError)
            {
                printf("  - First error at: n = %d, k = %d, error = %ld\n", n, k, error);
                firstError = 1;
            }
        }
    }

    printf("  - Maximum error: %ld\n", error);

    printf("  - Duration for %d computations of (%d,%d)\n", nRepeat, n_max, n_max/2);

    time_start = clock();
    for (i = 0; i < nRepeat; ++i)
        binomial_coefficient(n_max, n_max/2, 0);
    time_end = clock();
    printf("    logfact implementation: %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    time_start = clock();
    for (i = 0; i < nRepeat; ++i)
        binomial_coefficient(n_max, n_max/2, 1);
    time_end = clock();
    printf("    exact implementation: %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
}

/*!
 * Test the exactness of the full resolution wavelet transform in harmonic space.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_transform_axisym_lm_wav_test(int B, int L, int J_min, int seed)
{
  s2let_parameters_t parameters = {};
  parameters.B = B;
  parameters.L = L;
  parameters.J_min = J_min;

  clock_t time_start, time_end;
  double *wav_lm, *scal_lm;

  // Allocate the wavelet kernels
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, &parameters);

  // Compute the wavelet kernels
  time_start = clock();
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, &parameters);
  time_end = clock();
  printf("  - Generate wavelets  : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  complex double *f_wav_lm, *f_scal_lm, *flm, *flm_rec;
  s2let_allocate_lm(&flm, L);
  s2let_allocate_lm(&flm_rec, L);

  // Generate a random spherical harmonic decomposition
  s2let_lm_random_flm(flm, L, 0, seed);

  // Allocate space for the wavelet scales (their harmonic coefficients)
  s2let_transform_axisym_lm_allocate_f_wav(&f_wav_lm, &f_scal_lm, &parameters);

  // Perform the wavelet transform through exact harmonic tiling
  time_start = clock();
  s2let_transform_axisym_lm_wav_analysis(f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, &parameters);
  time_end = clock();
  printf("  - Wavelet analysis   : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Reconstruct the initial harmonic coefficients from those of the wavelets
  time_start = clock();
  s2let_transform_axisym_lm_wav_synthesis(flm_rec, f_wav_lm, f_scal_lm, wav_lm, scal_lm, &parameters);
  time_end = clock();
  printf("  - Wavelet synthesis  : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Compute the maximum absolute error on the harmonic coefficients
    /////////////////////////////////////////////////////////////////////////
  //  printf("  - size of flm :  %d\n",sizeof(flm));
 //   printf("  - size of flm_rec :  %d\n",sizeof(flm_rec));
 //   printf("  - analysis, flm[1] :  %6.5e+%6.5ei\n", flm[1]);
 //   printf("  - synthesis, flm_rec[1] :  %6.5e+%6.5ei\n", flm_rec[1]);
   /////////////////////////////////////////////////////////////////////////
  printf("  - Maximum abs error  : %6.5e\n",
	 maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

  free(flm);
  free(flm_rec);
  free(f_wav_lm);
  free(f_scal_lm);
  free(wav_lm);
  free(scal_lm);
}

/*!
 * Test the exactness of the multiresolution wavelet transform in harmonic space.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_transform_axisym_lm_wav_multires_test(int B, int L, int J_min, int seed)
{
  s2let_parameters_t parameters = {};
  parameters.B = B;
  parameters.L = L;
  parameters.J_min = J_min;

  clock_t time_start, time_end;
  double *wav_lm, *scal_lm;

  // Allocate the wavelet kernels
  s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, &parameters);

  // Compute the wavelet kernels
  time_start = clock();
  s2let_transform_axisym_lm_wav(wav_lm, scal_lm, &parameters);
  time_end = clock();
  printf("  - Generate wavelets  : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  complex double *f_wav_lm, *f_scal_lm, *flm, *flm_rec;
  s2let_allocate_lm(&flm, L);
  s2let_allocate_lm(&flm_rec, L);

  // Generate a random spherical harmonic decomposition
  s2let_lm_random_flm(flm, L, 0, seed);

  // Allocate space for the wavelet scales (their harmonic coefficients)
  s2let_transform_axisym_lm_allocate_f_wav_multires(&f_wav_lm, &f_scal_lm, &parameters);

  // Perform the wavelet transform through exact harmonic tiling
  time_start = clock();
  s2let_transform_axisym_lm_wav_analysis_multires(f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, &parameters);
  time_end = clock();
  printf("  - Wavelet analysis   : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Reconstruct the initial harmonic coefficients from those of the wavelets
  time_start = clock();
  s2let_transform_axisym_lm_wav_synthesis_multires(flm_rec, f_wav_lm, f_scal_lm, wav_lm, scal_lm, &parameters);
  time_end = clock();
  printf("  - Wavelet synthesis  : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Compute the maximum absolute error on the harmonic coefficients
    /////////////////////////////////////////////////////////////////////////
   // printf("  - analysis, flm[1] :  %6.5e+%6.5ei\n", flm[1]);
   // printf("  - synthesis, flm_rec[1] :  %6.5e+%6.5ei\n", flm_rec[1]);
    /////////////////////////////////////////////////////////////////////////
  printf("  - Maximum abs error  : %6.5e\n",
	 maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

  free(flm);
  free(flm_rec);
  free(f_wav_lm);
  free(f_scal_lm);
  free(wav_lm);
  free(scal_lm);
}

/*!
 * Test the exactness of the full resolution directional wavelet transform
 * in harmonic space.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_wav_transform_harmonic_test(int B, int L, int J_min, int N, int spin, int seed)
{
  s2let_parameters_t parameters = {};
  parameters.B = B;
  parameters.L = L;
  parameters.J_min = J_min;
  parameters.N = N;
  parameters.spin = spin;
  parameters.upsample = 1;
  parameters.normalization = S2LET_WAV_NORM_DEFAULT;
  parameters.original_spin = 0;

  clock_t time_start, time_end;
  complex double *psi;
  double *phi;

  // Allocate the wavelet kernels
  s2let_tiling_wavelet_allocate(&psi, &phi, &parameters);

  // Compute the wavelet kernels
  time_start = clock();
  s2let_tiling_wavelet(psi, phi, &parameters);
  time_end = clock();
  printf("  - Generate wavelets  : %4.4f seconds\n",
     (time_end - time_start) / (double)CLOCKS_PER_SEC);

  complex double *f_wav_lmn, *f_scal_lm, *flm, *flm_rec;
  s2let_allocate_lm(&flm, L);
  s2let_allocate_lm(&flm_rec, L);

  // Generate a random spherical harmonic decomposition
  s2let_lm_random_flm(flm, L, spin, seed);
    /////////////////////////////////////////////////////////////////////////////////////
    // printf("  -Generated random flm [1]:  %e+%ei\n",flm[1]);
    /////////////////////////////////////////////////////////////////////////////////////

  // Allocate space for the wavelet scales (their harmonic/Wigner coefficients)
  s2let_allocate_lmn_f_wav(&f_wav_lmn, &f_scal_lm, &parameters);

  // Perform the wavelet transform through exact harmonic tiling
  time_start = clock();
  s2let_analysis_lm2lmn(f_wav_lmn, f_scal_lm, flm, psi, phi, &parameters);
  time_end = clock();
  printf("  - Wavelet analysis   : %4.4f seconds\n",
     (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Reconstruct the initial harmonic coefficients from those of the wavelets
  time_start = clock();
  s2let_synthesis_lmn2lm(flm_rec, f_wav_lmn, f_scal_lm, psi, phi, &parameters);
  time_end = clock();
  printf("  - Wavelet synthesis  : %4.4f seconds\n",
     (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Compute the maximum absolute error on the harmonic coefficients
   /////////////////////////////////////////////////////////////////////////////////////
   // printf("  - size of flm :  %d\n",sizeof(flm));
   // printf("  - size of flm_rec :  %d\n",sizeof(flm_rec));
   // printf("  - analysis, flm[8] :  %6.5e+%6.5ei\n",flm[8]);
   // printf("  - synthesis, flm_rec[8] :  %6.5e+%6.5ei\n",flm_rec[8]);
    /////////////////////////////////////////////////////////////////////////////////////
  printf("  - Maximum abs error  : %6.5e\n",
     maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

  free(flm);
  free(flm_rec);
  free(f_wav_lmn);
  free(f_scal_lm);
  free(psi);
  free(phi);
}

/*!
 * Test the exactness of the multi-resolution directional wavelet transform
 * in harmonic space.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_wav_transform_harmonic_multires_test(int B, int L, int J_min, int N, int spin, int seed)
{
  s2let_parameters_t parameters = {};
  parameters.B = B;
  parameters.L = L;
  parameters.J_min = J_min;
  parameters.N = N;
  parameters.spin = spin;
  parameters.normalization = S2LET_WAV_NORM_DEFAULT;
  parameters.original_spin = 0;

  clock_t time_start, time_end;
  complex double *psi;
  double *phi;

  // Allocate the wavelet kernels
  s2let_tiling_wavelet_allocate(&psi, &phi, &parameters);

  // Compute the wavelet kernels
  time_start = clock();
  s2let_tiling_wavelet(psi, phi, &parameters);
  time_end = clock();
  printf("  - Generate wavelets  : %4.4f seconds\n",
     (time_end - time_start) / (double)CLOCKS_PER_SEC);

  complex double *f_wav_lmn, *f_scal_lm, *flm, *flm_rec;
  s2let_allocate_lm(&flm, L);
  s2let_allocate_lm(&flm_rec, L);

  // Generate a random spherical harmonic decomposition
  s2let_lm_random_flm(flm, L, spin, seed);

  // Allocate space for the wavelet scales (their harmonic/Wigner coefficients)
  s2let_allocate_lmn_f_wav(&f_wav_lmn, &f_scal_lm, &parameters);
/////////////////////////////////////////////////////////////////////////////////////
    // printf("  -analysis, f_wav_lmn :  %e+%ei\n",f_wav_lmn);
/////////////////////////////////////////////////////////////////////////////////////

  // Perform the wavelet transform through exact harmonic tiling
  time_start = clock();
  s2let_analysis_lm2lmn(f_wav_lmn, f_scal_lm, flm, psi, phi, &parameters);
  time_end = clock();
  printf("  - Wavelet analysis   : %4.4f seconds\n",
     (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Reconstruct the initial harmonic coefficients from those of the wavelets
  time_start = clock();
  s2let_synthesis_lmn2lm(flm_rec, f_wav_lmn, f_scal_lm, psi, phi, &parameters);
  time_end = clock();
  printf("  - Wavelet synthesis  : %4.4f seconds\n",
     (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Compute the maximum absolute error on the harmonic coefficients
    /////////////////////////////////////////////////////////////////////////////////////
   // printf("  -analysis, flm[1] :  %6.5e+%6.5ei\n",flm[1]);
   // printf("  -synthesis, flm_rec[1] :  %6.5e+%6.5ei\n",flm_rec[1]);
    /////////////////////////////////////////////////////////////////////////////////////
  printf("  - Maximum abs error  : %6.5e\n",
     maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

  free(flm);
  free(flm_rec);
  free(f_wav_lmn);
  free(f_scal_lm);
  free(psi);
  free(phi);
}

// CUREVELETS: 
/*!
 * Test the exactness of the full resolution curvelet transform
 * in harmonic space.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_cur_transform_harmonic_test(int B, int L, int J_min, int N, int spin, int seed)
{
    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = spin;
    parameters.upsample = 1;
    parameters.normalization = S2LET_WAV_NORM_DEFAULT;
    parameters.original_spin = 0;
    
    clock_t time_start, time_end;
    complex double *psi;
    double *phi;
    
    // Allocate the curvelet kernels
    s2let_tiling_curvelet_allocate(&psi, &phi, &parameters);
    
    // Compute the curvelet kernels
    time_start = clock();
    s2let_tiling_curvelet(psi, phi, &parameters);
    time_end = clock();
    printf("  - Generate curvelets  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    complex double *f_cur_lmn, *f_scal_lm, *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    
    // Generate a random spherical harmonic decomposition
    s2let_lm_random_flm(flm, L, spin, seed);
    /////////////////////////////////////////////////////////////////////////////////////
    // printf("  -Generated random flm :  %e+%ei\n",flm);
    /////////////////////////////////////////////////////////////////////////////////////
    
    // Allocate space for the curvelet scales (their harmonic/Wigner coefficients)
    s2let_allocate_lmn_f_cur(&f_cur_lmn, &f_scal_lm, &parameters);
    
    // Perform the curvelet transform through exact harmonic tiling
    time_start = clock();
    s2let_analysis_cur_lm2lmn(f_cur_lmn, f_scal_lm, flm, psi, phi, &parameters);
    time_end = clock();
    /////////////////////////////////////////////////////////////////////////////////////
    // printf("  -analysis, f_cur_lmn :  %e+%ei\n",f_cur_lmn);
    /////////////////////////////////////////////////////////////////////////////////////
    time_end = clock();
    printf("  - Curvelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Reconstruct the initial harmonic coefficients from those of the curvelets
    time_start = clock();
    s2let_synthesis_cur_lmn2lm(flm_rec, f_cur_lmn, f_scal_lm, psi, phi, &parameters);
    time_end = clock();
/////////////////////////////////////////////////////////////////////////////////////
   // printf("  - size of flm :  %d\n",sizeof(flm));
   // printf("  - size of flm_rec :  %d\n",sizeof(flm_rec));
   // printf("  - analysis, flm[1] :  %6.5e+%ei\n", flm[1]);
   // printf("  - synthesis, flm_rec[1] :  %6.5e+%ei\n", flm_rec[1]);
/////////////////////////////////////////////////////////////////////////////////////
    printf("  - Curvelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);
    
    free(flm);
    free(flm_rec);
    free(f_cur_lmn);
    free(f_scal_lm);
    free(psi);
    free(phi);
}

/*!
 * Test the exactness of the multi-resolution curvelet transform
 * in harmonic space.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_cur_transform_harmonic_multires_test(int B, int L, int J_min, int N, int spin, int seed)
{
    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = spin;
    parameters.normalization = S2LET_WAV_NORM_DEFAULT;
    parameters.original_spin = 0;
    
    clock_t time_start, time_end;
    complex double *psi;
    double *phi;
    
    // Allocate the wavelet kernels
    s2let_tiling_curvelet_allocate(&psi, &phi, &parameters);
    
    // Compute the wavelet kernels
    time_start = clock();
    s2let_tiling_curvelet(psi, phi, &parameters);
    time_end = clock();
    printf("  - Generate wavelets  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    complex double *f_cur_lmn, *f_scal_lm, *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    
    // Generate a random spherical harmonic decomposition
    s2let_lm_random_flm(flm, L, spin, seed);
    
    // Allocate space for the curvelet scales (their harmonic/Wigner coefficients)
    s2let_allocate_lmn_f_cur(&f_cur_lmn, &f_scal_lm, &parameters);
    
    // Perform the curvelet transform through exact harmonic tiling
    time_start = clock();
    s2let_analysis_cur_lm2lmn(f_cur_lmn, f_scal_lm, flm, psi, phi, &parameters);
    time_end = clock();
    printf("  - Curvelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Reconstruct the initial harmonic coefficients from those of the curvelets
    time_start = clock();
    s2let_synthesis_cur_lmn2lm(flm_rec, f_cur_lmn, f_scal_lm, psi, phi, &parameters);
    time_end = clock();
    printf("  - Curvelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Compute the maximum absolute error on the harmonic coefficients
    /////////////////////////////////////////////////////////////////////////////////////
    //printf("  - analysis, flm[1] :  %e+%ei\n", flm[1]);
    //printf("  - synthesis, flm_rec[1] :  %e+%ei\n", flm_rec[1]);
    /////////////////////////////////////////////////////////////////////////////////////
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);
    
    free(flm);
    free(flm_rec);
    free(f_cur_lmn);
    free(f_scal_lm);
    free(psi);
    free(phi);
}




/*!
 * Test the exactness of the full resolution wavelet transform in real space for complex functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_transform_axisym_wav_test(int B, int L, int J_min, int seed)
{
  s2let_parameters_t parameters = {};
  parameters.B = B;
  parameters.L = L;
  parameters.J_min = J_min;

  clock_t time_start, time_end;
  int spin = 0;
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;
  //int J = s2let_j_max(L, B);

  complex double *f, *f_rec, *flm, *flm_rec;
  s2let_allocate_lm(&flm, L);
  s2let_allocate_lm(&flm_rec, L);
  s2let_allocate_mw(&f, L);
  s2let_allocate_mw(&f_rec, L);

  // Generate random harmonic coefficients for a complex signal
  s2let_lm_random_flm(flm, L, 0, seed);

  // Construct the corresponding signal on the sphere (MW sampling)
  ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

  // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
  complex double *f_wav, *f_scal;
  s2let_transform_axisym_allocate_mw_f_wav(&f_wav, &f_scal, &parameters);

  // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
  time_start = clock();
  s2let_transform_axisym_wav_analysis_mw(f_wav, f_scal, f, &parameters);
  time_end = clock();
  printf("  - Wavelet analysis   : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Reconstruct the initial signal from the wavelet maps from scratch
  time_start = clock();
  s2let_transform_axisym_wav_synthesis_mw(f_rec, f_wav, f_scal, &parameters);
  time_end = clock();
  printf("  - Wavelet synthesis  : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Compute the initial harmonic coefficients back
  ssht_core_mw_forward_sov_conv_sym(flm_rec, f_rec, L, spin, dl_method, verbosity);

  // Compute the maximum absolute error on the harmonic coefficients
  printf("  - Maximum abs error  : %6.5e\n",
	 maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

  free(f);
  free(f_rec);
  free(f_wav);
  free(f_scal);
}

/*!
 * Test the exactness of the full resolution wavelet transform in real space for real functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_transform_axisym_wav_real_test(int B, int L, int J_min, int seed)
{
  s2let_parameters_t parameters = {};
  parameters.B = B;
  parameters.L = L;
  parameters.J_min = J_min;

  clock_t time_start, time_end;
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;
  //int J = s2let_j_max(L, B);

  complex *flm, *flm_rec;
  double *f, *f_rec;
  s2let_allocate_lm(&flm, L);
  s2let_allocate_lm(&flm_rec, L);
  s2let_allocate_mw_real(&f, L);
  s2let_allocate_mw_real(&f_rec, L);

  // Generate random harmonic coefficients for a real signal
  s2let_lm_random_flm_real(flm, L, seed);

  // Construct the corresponding signal on the sphere (MW sampling)
  ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);

  // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
  double *f_wav, *f_scal;
  s2let_transform_axisym_allocate_mw_f_wav_real(&f_wav, &f_scal, &parameters);

  // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
  time_start = clock();
  s2let_transform_axisym_wav_analysis_mw_real(f_wav, f_scal, f, &parameters);
  time_end = clock();
  printf("  - Wavelet analysis   : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Reconstruct the initial signal from the wavelet maps from scratch
  time_start = clock();
  s2let_transform_axisym_wav_synthesis_mw_real(f_rec, f_wav, f_scal, &parameters);
  time_end = clock();
  printf("  - Wavelet synthesis  : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Compute the initial harmonic coefficients back
  ssht_core_mw_forward_sov_conv_sym_real(flm_rec, f_rec, L, dl_method, verbosity);

  // Compute the maximum absolute error on the harmonic coefficients
  printf("  - Maximum abs error  : %6.5e\n",
	 maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

  free(f);
  free(f_rec);
  free(f_wav);
  free(f_scal);
}


/*!
 * Test the exactness of the multiresolution wavelet transform in real space for complex functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_transform_axisym_wav_multires_test(int B, int L, int J_min, int seed)
{
  s2let_parameters_t parameters = {};
  parameters.B = B;
  parameters.L = L;
  parameters.J_min = J_min;

  clock_t time_start, time_end;
  int spin = 0;
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;
  //int J = s2let_j_max(L, B);

  complex double *f, *f_rec, *flm, *flm_rec;
  s2let_allocate_lm(&flm, L);
  s2let_allocate_lm(&flm_rec, L);
  s2let_allocate_mw(&f, L);
  s2let_allocate_mw(&f_rec, L);

  // Generate random harmonic coefficients for a complex signal
  s2let_lm_random_flm(flm, L, 0, seed);

  // Construct the corresponding signal on the sphere (MW sampling)
  ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

  // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
  complex double *f_wav, *f_scal;
  s2let_transform_axisym_allocate_mw_f_wav_multires(&f_wav, &f_scal, &parameters);

  // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
  time_start = clock();
  s2let_transform_axisym_wav_analysis_mw_multires(f_wav, f_scal, f, &parameters);
  time_end = clock();
  printf("  - Wavelet analysis   : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Reconstruct the initial signal from the wavelet maps from scratch
  time_start = clock();
  s2let_transform_axisym_wav_synthesis_mw_multires(f_rec, f_wav, f_scal, &parameters);
  time_end = clock();
  printf("  - Wavelet synthesis  : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Compute the initial harmonic coefficients back
  ssht_core_mw_forward_sov_conv_sym(flm_rec, f_rec, L, spin, dl_method, verbosity);

  // Compute the maximum absolute error on the harmonic coefficients
  printf("  - Maximum abs error  : %6.5e\n",
	 maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

  free(f);
  free(f_rec);
  free(f_wav);
  free(f_scal);
}


/*!
 * Test the exactness of the multiresolution wavelet transform in real space for real functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_transform_axisym_wav_multires_real_test(int B, int L, int J_min, int seed)
{
  s2let_parameters_t parameters = {};
  parameters.B = B;
  parameters.L = L;
  parameters.J_min = J_min;

  clock_t time_start, time_end;
  int verbosity = 0;
  ssht_dl_method_t dl_method = SSHT_DL_RISBO;
  //int J = s2let_j_max(L, B);

  complex *flm, *flm_rec;
  double *f, *f_rec;
  s2let_allocate_lm(&flm, L);
  s2let_allocate_lm(&flm_rec, L);
  s2let_allocate_mw_real(&f, L);
  s2let_allocate_mw_real(&f_rec, L);

  // Generate random harmonic coefficients for a real signal
  s2let_lm_random_flm_real(flm, L, seed);

  // Construct the corresponding signal on the sphere (MW sampling)
  ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);

  // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
  double *f_wav, *f_scal;
  s2let_transform_axisym_allocate_mw_f_wav_multires_real(&f_wav, &f_scal, &parameters);

  // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
  time_start = clock();
  s2let_transform_axisym_wav_analysis_mw_multires_real(f_wav, f_scal, f, &parameters);
  time_end = clock();
  printf("  - Wavelet analysis   : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Reconstruct the initial signal from the wavelet maps from scratch
  time_start = clock();
  s2let_transform_axisym_wav_synthesis_mw_multires_real(f_rec, f_wav, f_scal, &parameters);
  time_end = clock();
  printf("  - Wavelet synthesis  : %4.4f seconds\n",
	 (time_end - time_start) / (double)CLOCKS_PER_SEC);

  // Compute the initial harmonic coefficients back
  ssht_core_mw_forward_sov_conv_sym_real(flm_rec, f_rec, L, dl_method, verbosity);

  // Compute the maximum absolute error on the harmonic coefficients
  printf("  - Maximum abs error  : %6.5e\n",
	 maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

  free(f);
  free(f_rec);
  free(f_wav);
  free(f_scal);
}

/*!
 * Test the exactness of the full resolution directional wavelet transform
 * in pixel space for complex functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_wav_transform_mw_test(int B, int L, int J_min, int N, int spin, int seed)
{
    clock_t time_start, time_end;

    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = spin;
    parameters.upsample = 1;
    parameters.normalization = S2LET_WAV_NORM_DEFAULT;
    parameters.original_spin = 0;

    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);

    complex double *f, *f_rec, *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    s2let_allocate_mw(&f, L);
    s2let_allocate_mw(&f_rec, L);

    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm(flm, L, spin, seed);
    // For debugging:
    int arrayind, arrayind_min;
    FILE *fp9;
    fp9=fopen("1_wav_flm_randgen_mw_test.dat", "w");
    arrayind_min= spin*spin;
    for (arrayind = arrayind_min; arrayind < L*L; arrayind++ )
    {
        fprintf(fp9, "%f, %f\n", creal(flm[arrayind]), cimag(flm[arrayind]));
    }
    fclose(fp9);
    
    
    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);
    // For debugging:
    // Open data file '"f_rec_afterssht_mw_inverse.dat"' to write out f
    FILE *fp2;
    arrayind =0 ;
    fp2=fopen("2_wav_f_afterssht_mw_inverse.dat", "w");
    for (arrayind =arrayind_min; arrayind < L*L; arrayind++ ) //sizeof(f)*2
    {
        fprintf(fp2, "%f, %f\n", creal(f[arrayind]), cimag(f[arrayind]));
    }
     fclose(fp2);
    

    // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
    complex double *f_wav, *f_scal;
    s2let_allocate_f_wav(&f_wav, &f_scal, &parameters);

    // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_px2wav(f_wav, f_scal, f, &parameters);
    time_end = clock();
    printf("  - Wavelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    // Open data file '"f_cur_afters2let_ana_px2cur.dat"' to write out f_cur
    FILE *fp3;
    arrayind =0 ;
    fp3=fopen("3b_f_wav_afters2let_ana_px2wav.dat", "w");
    for (arrayind =arrayind_min; arrayind < sizeof(f_wav)*2; arrayind++ )
    {
        fprintf(fp3, "%f, %f\n", creal(f_wav[arrayind]), cimag(f_wav[arrayind]));
    }
         fclose(fp3);
    
    
    
    // Reconstruct the initial signal from the wavelet maps from scratch
    time_start = clock();
    s2let_synthesis_wav2px(f_rec, f_wav, f_scal, &parameters);
    time_end = clock();
    printf("  - Wavelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    // For debugging:
    // Open data file '"f_cur_afters2let_ana_px2cur.dat"' to write out f_cur
    FILE *fp4;
    arrayind =0 ;
    fp4=fopen("5_wav_f_rec_afters2let_syn_wav2px.dat", "w");
    for (arrayind =arrayind_min; arrayind < L*L; arrayind++ ) //sizeof(f_rec)*2
    {
        fprintf(fp4, "%f, %f\n", creal(f_rec[arrayind]), cimag(f_rec[arrayind]));
    }
     fclose(fp4);
    
    
    
    
    // Convert back to harmonic coefficients
    ssht_core_mw_forward_sov_conv_sym(flm_rec, f_rec, L, spin, dl_method, verbosity);
    // For debugging:
    // Open data file '"flm_rec_afterssht_mw_forward.dat"' to write out flm_rec
    FILE *fp10;
    arrayind =0 ;
    fp10=fopen("6_wav_flm_rec_afterssht_mw_forward.dat", "w");
    for (arrayind =arrayind_min; arrayind < L*L; arrayind++ )  //sizeof(flm_rec)*2
    {
       fprintf(fp10, "%f, %f\n", creal(flm_rec[arrayind]), cimag(flm_rec[arrayind]));
    }
    fclose(fp10);

    
    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

    free(f);
    free(f_rec);
    free(flm);
    free(flm_rec);
    free(f_wav);
    free(f_scal);
}

/*!
 * Test the exactness of the full resolution directional wavelet transform
 * in pixel space for real functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_wav_transform_mw_real_test(int B, int L, int J_min, int N, int seed)
{
    clock_t time_start, time_end;

    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = 0;
    parameters.upsample = 1;

    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);

    double *f, *f_rec;
    complex double *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    s2let_allocate_mw_real(&f, L);
    s2let_allocate_mw_real(&f_rec, L);

    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm_real(flm, L, seed);

    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);

    // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
    double *f_wav, *f_scal;
    s2let_allocate_f_wav_real(&f_wav, &f_scal, &parameters);

    // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_px2wav_real(f_wav, f_scal, f, &parameters);
    time_end = clock();
    printf("  - Wavelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Reconstruct the initial signal from the wavelet maps from scratch
    time_start = clock();
    s2let_synthesis_wav2px_real(f_rec, f_wav, f_scal, &parameters);
    time_end = clock();
    printf("  - Wavelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Convert back to harmonic coefficients
    ssht_core_mw_forward_sov_conv_sym_real(flm_rec, f_rec, L, dl_method, verbosity);

    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

    free(f);
    free(f_rec);
    free(flm);
    free(flm_rec);
    free(f_wav);
    free(f_scal);
}

/*!
 * Test the exactness of the multi-resolution directional wavelet transform
 * in pixel space for complex functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_wav_transform_mw_multires_test(int B, int L, int J_min, int N, int spin, int seed)
{
    clock_t time_start, time_end;

    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = spin;
    parameters.normalization = S2LET_WAV_NORM_DEFAULT;
    parameters.original_spin = 0;

    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);

    complex double *f, *f_rec, *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    s2let_allocate_mw(&f, L);
    s2let_allocate_mw(&f_rec, L);

    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm(flm, L, spin, seed);

    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

    // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
    complex double *f_wav, *f_scal;
    s2let_allocate_f_wav(&f_wav, &f_scal, &parameters);

    // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_px2wav(f_wav, f_scal, f, &parameters);
    time_end = clock();
    printf("  - Wavelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Reconstruct the initial signal from the wavelet maps from scratch
    time_start = clock();
    s2let_synthesis_wav2px(f_rec, f_wav, f_scal, &parameters);
    time_end = clock();
    printf("  - Wavelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Convert back to harmonic coefficients
    ssht_core_mw_forward_sov_conv_sym(flm_rec, f_rec, L, spin, dl_method, verbosity);

    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

    free(f);
    free(f_rec);
    free(flm);
    free(flm_rec);
    free(f_wav);
    free(f_scal);
}

/*!
 * Test the exactness of the multi-resolution directional wavelet transform
 * in pixel space for real functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_wav_transform_mw_multires_real_test(int B, int L, int J_min, int N, int seed)
{
    clock_t time_start, time_end;

    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = 0;
    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);

    double *f, *f_rec;
    complex double *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    s2let_allocate_mw_real(&f, L);
    s2let_allocate_mw_real(&f_rec, L);

    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm_real(flm, L, seed);

    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);

    // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
    double *f_wav, *f_scal;
    s2let_allocate_f_wav_real(&f_wav, &f_scal, &parameters);

    // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_px2wav_real(f_wav, f_scal, f, &parameters);
    time_end = clock();
    printf("  - Wavelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Reconstruct the initial signal from the wavelet maps from scratch
    time_start = clock();
    s2let_synthesis_wav2px_real(f_rec, f_wav, f_scal, &parameters);
    time_end = clock();
    printf("  - Wavelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Convert back to harmonic coefficients
    ssht_core_mw_forward_sov_conv_sym_real(flm_rec, f_rec, L, dl_method, verbosity);

    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

    free(f);
    free(f_rec);
    free(flm);
    free(flm_rec);
    free(f_wav);
    free(f_scal);
}

/*!
 * Test the exactness of the full resolution directional wavelet transform
 * in pixel space for complex functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_wav_transform_mwss_test(int B, int L, int J_min, int N, int spin, int seed)
{
    clock_t time_start, time_end;

    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = spin;
    parameters.upsample = 1;
    parameters.normalization = S2LET_WAV_NORM_DEFAULT;
    parameters.original_spin = 0;
    parameters.sampling_scheme = S2LET_SAMPLING_MW_SS;

    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);

    complex double *f, *f_rec, *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    s2let_allocate_mwss(&f, L);
    s2let_allocate_mwss(&f_rec, L);

    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm(flm, L, spin, seed);

    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym_ss(f, flm, L, spin, dl_method, verbosity);

    // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
    complex double *f_wav, *f_scal;
    s2let_allocate_f_wav(&f_wav, &f_scal, &parameters);

    // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_px2wav(f_wav, f_scal, f, &parameters);
    time_end = clock();
    printf("  - Wavelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Reconstruct the initial signal from the wavelet maps from scratch
    time_start = clock();
    s2let_synthesis_wav2px(f_rec, f_wav, f_scal, &parameters);
    time_end = clock();
    printf("  - Wavelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Convert back to harmonic coefficients
    ssht_core_mw_forward_sov_conv_sym_ss(flm_rec, f_rec, L, spin, dl_method, verbosity);

    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

    free(f);
    free(f_rec);
    free(flm);
    free(flm_rec);
    free(f_wav);
    free(f_scal);
}

/*!
 * Test the exactness of the full resolution directional wavelet transform
 * in pixel space for real functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_wav_transform_mwss_real_test(int B, int L, int J_min, int N, int seed)
{
    clock_t time_start, time_end;

    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = 0;
    parameters.upsample = 1;
    parameters.sampling_scheme = S2LET_SAMPLING_MW_SS;

    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);

    double *f, *f_rec;
    complex double *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    s2let_allocate_mwss_real(&f, L);
    s2let_allocate_mwss_real(&f_rec, L);

    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm_real(flm, L, seed);

    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym_ss_real(f, flm, L, dl_method, verbosity);

    // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
    double *f_wav, *f_scal;
    s2let_allocate_f_wav_real(&f_wav, &f_scal, &parameters);

    // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_px2wav_real(f_wav, f_scal, f, &parameters);
    time_end = clock();
    printf("  - Wavelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Reconstruct the initial signal from the wavelet maps from scratch
    time_start = clock();
    s2let_synthesis_wav2px_real(f_rec, f_wav, f_scal, &parameters);
    time_end = clock();
    printf("  - Wavelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Convert back to harmonic coefficients
    ssht_core_mw_forward_sov_conv_sym_ss_real(flm_rec, f_rec, L, dl_method, verbosity);

    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

    free(f);
    free(f_rec);
    free(flm);
    free(flm_rec);
    free(f_wav);
    free(f_scal);
}

/*!
 * Test the exactness of the multi-resolution directional wavelet transform
 * in pixel space for complex functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_wav_transform_mwss_multires_test(int B, int L, int J_min, int N, int spin, int seed)
{
    clock_t time_start, time_end;

    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = spin;
    parameters.normalization = S2LET_WAV_NORM_DEFAULT;
    parameters.original_spin = 0;
    parameters.sampling_scheme = S2LET_SAMPLING_MW_SS;

    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);

    complex double *f, *f_rec, *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    s2let_allocate_mwss(&f, L);
    s2let_allocate_mwss(&f_rec, L);

    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm(flm, L, spin, seed);

    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym_ss(f, flm, L, spin, dl_method, verbosity);

    // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
    complex double *f_wav, *f_scal;
    s2let_allocate_f_wav(&f_wav, &f_scal, &parameters);

    // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_px2wav(f_wav, f_scal, f, &parameters);
    time_end = clock();
    printf("  - Wavelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Reconstruct the initial signal from the wavelet maps from scratch
    time_start = clock();
    s2let_synthesis_wav2px(f_rec, f_wav, f_scal, &parameters);
    time_end = clock();
    printf("  - Wavelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Convert back to harmonic coefficients
    ssht_core_mw_forward_sov_conv_sym_ss(flm_rec, f_rec, L, spin, dl_method, verbosity);

    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

    free(f);
    free(f_rec);
    free(flm);
    free(flm_rec);
    free(f_wav);
    free(f_scal);
}

/*!
 * Test the exactness of the multi-resolution directional wavelet transform
 * in pixel space for real functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_wav_transform_mwss_multires_real_test(int B, int L, int J_min, int N, int seed)
{
    clock_t time_start, time_end;

    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = 0;
    parameters.sampling_scheme = S2LET_SAMPLING_MW_SS;

    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);

    double *f, *f_rec;
    complex double *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    s2let_allocate_mwss_real(&f, L);
    s2let_allocate_mwss_real(&f_rec, L);

    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm_real(flm, L, seed);

    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym_ss_real(f, flm, L, dl_method, verbosity);

    // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
    double *f_wav, *f_scal;
    s2let_allocate_f_wav_real(&f_wav, &f_scal, &parameters);

    // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_px2wav_real(f_wav, f_scal, f, &parameters);
    time_end = clock();
    printf("  - Wavelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Reconstruct the initial signal from the wavelet maps from scratch
    time_start = clock();
    s2let_synthesis_wav2px_real(f_rec, f_wav, f_scal, &parameters);
    time_end = clock();
    printf("  - Wavelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Convert back to harmonic coefficients
    ssht_core_mw_forward_sov_conv_sym_ss_real(flm_rec, f_rec, L, dl_method, verbosity);

    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

    free(f);
    free(f_rec);
    free(flm);
    free(flm_rec);
    free(f_wav);
    free(f_scal);
}

/*!
 * Test the exactness of the full resolution directional harmonic-to-wavelet
 * transform in pixel space for complex functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_wav_transform_lm2wav_test(int B, int L, int J_min, int N, int spin, int seed)
{
    clock_t time_start, time_end;

    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = spin;
    parameters.upsample = 1;
    parameters.normalization = S2LET_WAV_NORM_DEFAULT;
    parameters.original_spin = 0;
    parameters.verbosity = 0;
    parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);

    complex double *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);

    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm(flm, L, spin, seed);

    // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
    complex double *f_wav, *f_scal;
    s2let_allocate_f_wav(&f_wav, &f_scal, &parameters);

    // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_lm2wav(f_wav, f_scal, flm, &parameters);
    time_end = clock();
    printf("  - Wavelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Reconstruct the initial signal from the wavelet maps from scratch
    time_start = clock();
    s2let_synthesis_wav2lm(flm_rec, f_wav, f_scal, &parameters);
    time_end = clock();
    printf("  - Wavelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

    free(flm);
    free(flm_rec);
    free(f_wav);
    free(f_scal);
}

/*!
 * Test the exactness of the full resolution directional wavelet transform
 * in pixel space for real functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_wav_transform_lm2wav_real_test(int B, int L, int J_min, int N, int seed)
{
    clock_t time_start, time_end;

    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = 0;
    parameters.upsample = 1;

    parameters.verbosity = 0;
    parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);

    complex double *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);

    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm_real(flm, L, seed);

    // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
    double *f_wav, *f_scal;
    s2let_allocate_f_wav_real(&f_wav, &f_scal, &parameters);

    // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_lm2wav_real(f_wav, f_scal, flm, &parameters);
    time_end = clock();
    printf("  - Wavelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Reconstruct the initial signal from the wavelet maps from scratch
    time_start = clock();
    s2let_synthesis_wav2lm_real(flm_rec, f_wav, f_scal, &parameters);
    time_end = clock();
    printf("  - Wavelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

    free(flm);
    free(flm_rec);
    free(f_wav);
    free(f_scal);
}

/*!
 * Test the exactness of the multi-resolution directional wavelet transform
 * in pixel space for complex functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_wav_transform_lm2wav_multires_test(int B, int L, int J_min, int N, int spin, int seed)
{
    clock_t time_start, time_end;

    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = spin;
    parameters.normalization = S2LET_WAV_NORM_DEFAULT;
    parameters.original_spin = 0;
    parameters.verbosity = 0;
    parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);

    complex double *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);

    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm(flm, L, spin, seed);

    // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
    complex double *f_wav, *f_scal;
    s2let_allocate_f_wav(&f_wav, &f_scal, &parameters);

    // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_lm2wav(f_wav, f_scal, flm, &parameters);
    time_end = clock();
    printf("  - Wavelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Reconstruct the initial signal from the wavelet maps from scratch
    time_start = clock();
    s2let_synthesis_wav2lm(flm_rec, f_wav, f_scal, &parameters);
    time_end = clock();
    printf("  - Wavelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

    free(flm);
    free(flm_rec);
    free(f_wav);
    free(f_scal);
}

/*!
 * Test the exactness of the multi-resolution directional wavelet transform
 * in pixel space for real functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_wav_transform_lm2wav_multires_real_test(int B, int L, int J_min, int N, int seed)
{
    clock_t time_start, time_end;

    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = 0;
    parameters.verbosity = 0;
    parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);

    complex double *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);

    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm_real(flm, L, seed);

    // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
    double *f_wav, *f_scal;
    s2let_allocate_f_wav_real(&f_wav, &f_scal, &parameters);

    // Perform wavelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_lm2wav_real(f_wav, f_scal, flm, &parameters);
    time_end = clock();
    printf("  - Wavelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Reconstruct the initial signal from the wavelet maps from scratch
    time_start = clock();
    s2let_synthesis_wav2lm_real(flm_rec, f_wav, f_scal, &parameters);
    time_end = clock();
    printf("  - Wavelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);

    // Compute the maximum absolute error on the harmonic coefficients
    // printf("  - analysis, flm[1] :  %f+%fi\n", flm[1]);
    // printf("  - synthesis, flm_rec[1] : %f+%fi\n",flm_rec[1]);
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);

    free(flm);
    free(flm_rec);
    free(f_wav);
    free(f_scal);
}




// Curvelets:

/*!
 * Test the exactness of the full resolution curvelet transform
 * in pixel space for complex functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_cur_transform_mw_test(int B, int L, int J_min, int N, int spin, int seed)
{
    clock_t time_start, time_end;
    
    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = spin;
    parameters.upsample = 1;
    parameters.normalization = S2LET_WAV_NORM_DEFAULT;
    parameters.original_spin = 0;
    
    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);
    
    complex double *f, *f_rec, *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    s2let_allocate_mw(&f, L);
    s2let_allocate_mw(&f_rec, L);
    
    int arrayind;
    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm(flm, L, spin, seed);
    FILE *fp1;
    fp1=fopen("1_flm_randgen_mw_test.dat", "w");
    for (arrayind =0; arrayind < L*L; arrayind++ )   // sizeof(flm)*2
    {
        fprintf(fp1, "%f, %f\n", creal(flm[arrayind]), cimag(flm[arrayind]));
    }
    
    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);
    // For debugging:
    // Open data file '"f_rec_afterssht_mw_inverse.dat"' to write out f
    FILE *fp2;
    fp2=fopen("2_f_afterssht_mw_inverse.dat", "w");
    for (arrayind =0; arrayind <L*L ; arrayind++ )  //sizeof(f)*2
    {
        fprintf(fp2, "%f, %f\n", creal(f[arrayind]), cimag(f[arrayind]));
    }

    
    // Allocate space for curvelet maps on the sphere (corresponding to the triplet B/L/J_min)
    complex double *f_cur, *f_scal;
    s2let_allocate_f_cur(&f_cur, &f_scal, &parameters);
    
    // Perform curvelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_px2cur(f_cur, f_scal, f, &parameters);
    time_end = clock();
    printf("  - Curvelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    // For debugging:
    // Open data file '"f_cur_afters2let_ana_px2cur.dat"' to write out f_cur
    FILE *fp3;
    fp3=fopen("3b_f_cur_afters2let_ana_px2cur.dat", "w");
    for (arrayind =0; arrayind < L*L; arrayind++ ) //sizeof(f_cur)*2
    {
        fprintf(fp3, "%f, %f\n", creal(f_cur[arrayind]), cimag(f_cur[arrayind]));
    }
    
    
    // Reconstruct the initial signal from the curvelet maps from scratch
    time_start = clock();
    s2let_synthesis_cur2px(f_rec, f_cur, f_scal, &parameters);
    time_end = clock();
    printf("  - curvelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    // For debugging:
    // Open data file '"f_cur_afters2let_ana_px2cur.dat"' to write out f_cur
    FILE *fp4;
    fp4=fopen("5_f_rec_afters2let_syn_cur2px.dat", "w");
    for (arrayind =0; arrayind < L*L ; arrayind++ )//sizeof(f_rec)*2
    {
        fprintf(fp4, "%f, %f\n", creal(f_rec[arrayind]), cimag(f_rec[arrayind]));
    }
    
    
    // Convert back to harmonic coefficients
    ssht_core_mw_forward_sov_conv_sym(flm_rec, f_rec, L, spin, dl_method, verbosity);
    // For debugging:
    // Open data file '"flm_rec_afterssht_mw_forward.dat"' to write out flm_rec
    FILE *fp5;
    fp5=fopen("6_flm_rec_afterssht_mw_forward.dat", "w");
    for (arrayind =0; arrayind < L*L; arrayind++ ) //sizeof(flm_rec)*2
    {
        fprintf(fp5, "%f, %f\n", creal(flm_rec[arrayind]), cimag(flm_rec[arrayind]));
    }
    
    
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);
    fclose(fp5);
    
    
    // Compute the maximum absolute error on the harmonic coefficients
    /////////////////////////////////////////////////////////////////////////////////////
    printf("  - size of flm :  %ld\n",sizeof(flm));
    printf("  - size of flm_rec :  %ld\n",sizeof(flm_rec));
    printf("  - analysis, flm[1] :  %f, %f\n",creal(flm[1]), cimag(flm[1]));
    printf("  - synthesis, flm_rec[1] :   %f, %f\n",creal(flm_rec[1]), cimag(flm_rec[1]));
    printf("  - analysis, flm[5] :  %f, %f\n",creal(flm[5]), cimag(flm[5]));
    printf("  - synthesis, flm_rec[5] : %f, %f\n",creal(flm_rec[5]), cimag(flm_rec[5]));
  //  printf("  - analysis, flm[8] :  %6.5e+%6.5ei\n",flm[8]);
  //  printf("  - synthesis, flm_rec[8] :  %6.5e+%6.5ei\n",flm_rec[8]);
    /////////////////////////////////////////////////////////////////////////////////////
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);
    
    free(f);
    free(f_rec);
    free(flm);
    free(flm_rec);
    free(f_cur);
    free(f_scal);
}

/*!
 * Test the exactness of the full resolution curvelet transform
 * in pixel space for real functions.
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_cur_transform_mw_real_test(int B, int L, int J_min, int N, int seed)
{
    clock_t time_start, time_end;
    
    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = 0;
    parameters.upsample = 1;
    
    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);
    
    double *f, *f_rec;
    complex double *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    s2let_allocate_mw_real(&f, L);
    s2let_allocate_mw_real(&f_rec, L);
    
    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm_real(flm, L, seed);
    
    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);
    
    // Allocate space for curvelet maps on the sphere (corresponding to the triplet B/L/J_min)
    double *f_cur, *f_scal;
    s2let_allocate_f_cur_real(&f_cur, &f_scal, &parameters);
    
    // Perform curvelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_px2cur_real(f_cur, f_scal, f, &parameters);
    time_end = clock();
    printf("  - Curvelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Reconstruct the initial signal from the curvelet maps from scratch
    time_start = clock();
    s2let_synthesis_cur2px_real(f_rec, f_cur, f_scal, &parameters);
    time_end = clock();
    printf("  - Curvelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Convert back to harmonic coefficients
    ssht_core_mw_forward_sov_conv_sym_real(flm_rec, f_rec, L, dl_method, verbosity);
    
    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);
    
    free(f);
    free(f_rec);
    free(flm);
    free(flm_rec);
    free(f_cur);
    free(f_scal);
}

/*!
 * Test the exactness of the multi-resolution curvelet transform
 * in pixel space for complex functions.
 *
 * \param[in]  B curvelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First curvelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_cur_transform_mw_multires_test(int B, int L, int J_min, int N, int spin, int seed)
{
    clock_t time_start, time_end;
    
    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = spin;
    parameters.normalization = S2LET_WAV_NORM_DEFAULT;
    parameters.original_spin = 0;
    
    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);
    
    complex double *f, *f_rec, *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    s2let_allocate_mw(&f, L);
    s2let_allocate_mw(&f_rec, L);
    
    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm(flm, L, spin, seed);
    
    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);
    
    // Allocate space for curvelet maps on the sphere (corresponding to the triplet B/L/J_min)
    complex double *f_cur, *f_scal;
    s2let_allocate_f_cur(&f_cur, &f_scal, &parameters);
    
    // Perform curvelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_px2cur(f_cur, f_scal, f, &parameters);
    time_end = clock();
    printf("  - Curvelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Reconstruct the initial signal from the curvelet maps from scratch
    time_start = clock();
    s2let_synthesis_cur2px(f_rec, f_cur, f_scal, &parameters);
    time_end = clock();
    printf("  - Curvelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Convert back to harmonic coefficients
    ssht_core_mw_forward_sov_conv_sym(flm_rec, f_rec, L, spin, dl_method, verbosity);
    
    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);
    
    free(f);
    free(f_rec);
    free(flm);
    free(flm_rec);
    free(f_cur);
    free(f_scal);
}

/*!
 * Test the exactness of the multi-resolution curvelet transform
 * in pixel space for real functions.
 *
 * \param[in]  B curvelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First curvelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_cur_transform_mw_multires_real_test(int B, int L, int J_min, int N, int seed)
{
    clock_t time_start, time_end;
    
    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = 0;
    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);
    
    double *f, *f_rec;
    complex double *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    s2let_allocate_mw_real(&f, L);
    s2let_allocate_mw_real(&f_rec, L);
    
    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm_real(flm, L, seed);
    
    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);
    
    // Allocate space for curvelet maps on the sphere (corresponding to the triplet B/L/J_min)
    double *f_cur, *f_scal;
    s2let_allocate_f_cur_real(&f_cur, &f_scal, &parameters);
    
    // Perform curvelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_px2cur_real(f_cur, f_scal, f, &parameters);
    time_end = clock();
    printf("  - Curvelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Reconstruct the initial signal from the curvelet maps from scratch
    time_start = clock();
    s2let_synthesis_cur2px_real(f_rec, f_cur, f_scal, &parameters);
    time_end = clock();
    printf("  - Curvelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Convert back to harmonic coefficients
    ssht_core_mw_forward_sov_conv_sym_real(flm_rec, f_rec, L, dl_method, verbosity);
    
    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);
    
    free(f);
    free(f_rec);
    free(flm);
    free(flm_rec);
    free(f_cur);
    free(f_scal);
}

/*!
 * Test the exactness of the full resolution Curvelet transform
 * in pixel space for complex functions.
 *
 * \param[in]  B curvelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First curvelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_cur_transform_mwss_test(int B, int L, int J_min, int N, int spin, int seed)
{
    clock_t time_start, time_end;
    
    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = spin;
    parameters.upsample = 1;
    parameters.normalization = S2LET_WAV_NORM_DEFAULT;
    parameters.original_spin = 0;
    parameters.sampling_scheme = S2LET_SAMPLING_MW_SS;
    
    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);
    
    complex double *f, *f_rec, *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    s2let_allocate_mwss(&f, L);
    s2let_allocate_mwss(&f_rec, L);
    
    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm(flm, L, spin, seed);
    
    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym_ss(f, flm, L, spin, dl_method, verbosity);
    
    // Allocate space for curvelet maps on the sphere (corresponding to the triplet B/L/J_min)
    complex double *f_cur, *f_scal;
    s2let_allocate_f_cur(&f_cur, &f_scal, &parameters);
    
    // Perform curvelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_px2cur(f_cur, f_scal, f, &parameters);
    time_end = clock();
    printf("  - Curvelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Reconstruct the initial signal from the curvelet maps from scratch
    time_start = clock();
    s2let_synthesis_cur2px(f_rec, f_cur, f_scal, &parameters);
    time_end = clock();
    printf("  - Curvelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Convert back to harmonic coefficients
    ssht_core_mw_forward_sov_conv_sym_ss(flm_rec, f_rec, L, spin, dl_method, verbosity);
    
    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);
    
    free(f);
    free(f_rec);
    free(flm);
    free(flm_rec);
    free(f_cur);
    free(f_scal);
}

/*!
 * Test the exactness of the full resolution curvelet transform
 * in pixel space for real functions.
 *
 * \param[in]  B curvelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First curvelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_cur_transform_mwss_real_test(int B, int L, int J_min, int N, int seed)
{
    clock_t time_start, time_end;
    
    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = 0;
    parameters.upsample = 1;
    parameters.sampling_scheme = S2LET_SAMPLING_MW_SS;
    
    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);
    
    double *f, *f_rec;
    complex double *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    s2let_allocate_mwss_real(&f, L);
    s2let_allocate_mwss_real(&f_rec, L);
    
    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm_real(flm, L, seed);
    
    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym_ss_real(f, flm, L, dl_method, verbosity);
    
    // Allocate space for curvelet maps on the sphere (corresponding to the triplet B/L/J_min)
    double *f_cur, *f_scal;
    s2let_allocate_f_cur_real(&f_cur, &f_scal, &parameters);
    
    // Perform curvelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_px2cur_real(f_cur, f_scal, f, &parameters);
    time_end = clock();
    printf("  - Curvelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Reconstruct the initial signal from the curvelet maps from scratch
    time_start = clock();
    s2let_synthesis_cur2px_real(f_rec, f_cur, f_scal, &parameters);
    time_end = clock();
    printf("  - Curvelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Convert back to harmonic coefficients
    ssht_core_mw_forward_sov_conv_sym_ss_real(flm_rec, f_rec, L, dl_method, verbosity);
    
    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);
    
    free(f);
    free(f_rec);
    free(flm);
    free(flm_rec);
    free(f_cur);
    free(f_scal);
}

/*!
 * Test the exactness of the multi-resolution curvelet transform
 * in pixel space for complex functions.
 *
 * \param[in]  B curvelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First curvelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_cur_transform_mwss_multires_test(int B, int L, int J_min, int N, int spin, int seed)
{
    clock_t time_start, time_end;
    
    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = spin;
    parameters.normalization = S2LET_WAV_NORM_DEFAULT;
    parameters.original_spin = 0;
    parameters.sampling_scheme = S2LET_SAMPLING_MW_SS;
    
    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);
    
    complex double *f, *f_rec, *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    s2let_allocate_mwss(&f, L);
    s2let_allocate_mwss(&f_rec, L);
    
    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm(flm, L, spin, seed);
    
    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym_ss(f, flm, L, spin, dl_method, verbosity);
    
    // Allocate space for curvelet maps on the sphere (corresponding to the triplet B/L/J_min)
    complex double *f_cur, *f_scal;
    s2let_allocate_f_cur(&f_cur, &f_scal, &parameters);
    
    // Perform curvelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_px2cur(f_cur, f_scal, f, &parameters);
    time_end = clock();
    printf("  - Curvelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Reconstruct the initial signal from the curvelet maps from scratch
    time_start = clock();
    s2let_synthesis_cur2px(f_rec, f_cur, f_scal, &parameters);
    time_end = clock();
    printf("  - Curvelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Convert back to harmonic coefficients
    ssht_core_mw_forward_sov_conv_sym_ss(flm_rec, f_rec, L, spin, dl_method, verbosity);
    
    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);
    
    free(f);
    free(f_rec);
    free(flm);
    free(flm_rec);
    free(f_cur);
    free(f_scal);
}

/*!
 * Test the exactness of the multi-resolution curvelet transform
 * in pixel space for real functions.
 *
 * \param[in]  B curvelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First curvelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_cur_transform_mwss_multires_real_test(int B, int L, int J_min, int N, int seed)
{
    clock_t time_start, time_end;
    
    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = 0;
    parameters.sampling_scheme = S2LET_SAMPLING_MW_SS;
    
    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);
    
    double *f, *f_rec;
    complex double *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    s2let_allocate_mwss_real(&f, L);
    s2let_allocate_mwss_real(&f_rec, L);
    
    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm_real(flm, L, seed);
    
    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym_ss_real(f, flm, L, dl_method, verbosity);
    
    // Allocate space for curvelet maps on the sphere (corresponding to the triplet B/L/J_min)
    double *f_cur, *f_scal;
    s2let_allocate_f_cur_real(&f_cur, &f_scal, &parameters);
    
    // Perform curvelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_px2cur_real(f_cur, f_scal, f, &parameters);
    time_end = clock();
    printf("  - Curvelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Reconstruct the initial signal from the curvelet maps from scratch
    time_start = clock();
    s2let_synthesis_cur2px_real(f_rec, f_cur, f_scal, &parameters);
    time_end = clock();
    printf("  - Curvelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Convert back to harmonic coefficients
    ssht_core_mw_forward_sov_conv_sym_ss_real(flm_rec, f_rec, L, dl_method, verbosity);
    
    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);
    
    free(f);
    free(f_rec);
    free(flm);
    free(flm_rec);
    free(f_cur);
    free(f_scal);
}

/*!
 * Test the exactness of the full resolution harmonic-to-curvelet
 * transform in pixel space for complex functions.
 *
 * \param[in]  B curvelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First curvelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_cur_transform_lm2cur_test(int B, int L, int J_min, int N, int spin, int seed)
{
    clock_t time_start, time_end;
    
    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = spin;
    parameters.upsample = 1;
    parameters.normalization = S2LET_WAV_NORM_DEFAULT;
    parameters.original_spin = 0;
    parameters.verbosity = 0;
    parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);
    
    complex double *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    
    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm(flm, L, spin, seed);
    
    // Allocate space for curvelet maps on the sphere (corresponding to the triplet B/L/J_min)
    complex double *f_cur, *f_scal;
    s2let_allocate_f_cur(&f_cur, &f_scal, &parameters);
    
    // Perform curvelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_lm2cur(f_cur, f_scal, flm, &parameters);
    time_end = clock();
    printf("  - Curvelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Reconstruct the initial signal from the curvelet maps from scratch
    time_start = clock();
    s2let_synthesis_cur2lm(flm_rec, f_cur, f_scal, &parameters);
    time_end = clock();
    printf("  - Curvelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);
    
    free(flm);
    free(flm_rec);
    free(f_cur);
    free(f_scal);
}

/*!
 * Test the exactness of the full resolution curvelet transform
 * in pixel space for real functions.
 *
 * \param[in]  B curvelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First curvelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_cur_transform_lm2cur_real_test(int B, int L, int J_min, int N, int seed)
{
    clock_t time_start, time_end;
    
    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = 0;
    parameters.upsample = 1;
    
    parameters.verbosity = 0;
    parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);
    
    complex double *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    
    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm_real(flm, L, seed);
    
    // Allocate space for curvelet maps on the sphere (corresponding to the triplet B/L/J_min)
    double *f_cur, *f_scal;
    s2let_allocate_f_cur_real(&f_cur, &f_scal, &parameters);
    
    // Perform curvelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_lm2cur_real(f_cur, f_scal, flm, &parameters);
    time_end = clock();
    printf("  - curvelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Reconstruct the initial signal from the curvelet maps from scratch
    time_start = clock();
    s2let_synthesis_cur2lm_real(flm_rec, f_cur, f_scal, &parameters);
    time_end = clock();
    printf("  - curvelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);
    
    free(flm);
    free(flm_rec);
    free(f_cur);
    free(f_scal);
}

/*!
 * Test the exactness of the multi-resolution directional curvelet transform
 * in pixel space for complex functions.
 *
 * \param[in]  B curvelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First curvelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_cur_transform_lm2cur_multires_test(int B, int L, int J_min, int N, int spin, int seed)
{
    clock_t time_start, time_end;
    
    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = spin;
    parameters.normalization = S2LET_WAV_NORM_DEFAULT;
    parameters.original_spin = 0;
    parameters.verbosity = 0;
    parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);
    
    complex double *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    
    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm(flm, L, spin, seed);
    
    // Allocate space for curvelet maps on the sphere (corresponding to the triplet B/L/J_min)
    complex double *f_cur, *f_scal;
    s2let_allocate_f_cur(&f_cur, &f_scal, &parameters);
    
    // Perform curvelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_lm2cur(f_cur, f_scal, flm, &parameters);
    time_end = clock();
    printf("  - curvelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Reconstruct the initial signal from the curvelet maps from scratch
    time_start = clock();
    s2let_synthesis_cur2lm(flm_rec, f_cur, f_scal, &parameters);
    time_end = clock();
    printf("  - curvelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Compute the maximum absolute error on the harmonic coefficients
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);
    
    free(flm);
    free(flm_rec);
    free(f_cur);
    free(f_scal);
}

/*!
 * Test the exactness of the multi-resolution directional curvelet transform
 * in pixel space for real functions.
 *
 * \param[in]  B curvelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First curvelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \param[in]  spin Spin number.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_cur_transform_lm2cur_multires_real_test(int B, int L, int J_min, int N, int seed)
{
    clock_t time_start, time_end;
    
    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = N;
    parameters.spin = 0;
    parameters.verbosity = 0;
    parameters.dl_method = SSHT_DL_RISBO;
    //int J = s2let_j_max(L, B);
    
    complex double *flm, *flm_rec;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_lm(&flm_rec, L);
    
    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm_real(flm, L, seed);
    
    // Allocate space for curvelet maps on the sphere (corresponding to the triplet B/L/J_min)
    double *f_cur, *f_scal;
    s2let_allocate_f_cur_real(&f_cur, &f_scal, &parameters);
    
    // Perform curvelet analysis from scratch with all signals given on the sphere (MW sampling)
    time_start = clock();
    s2let_analysis_lm2cur_real(f_cur, f_scal, flm, &parameters);
    time_end = clock();
    printf("  - curvelet analysis   : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Reconstruct the initial signal from the curvelet maps from scratch
    time_start = clock();
    s2let_synthesis_cur2lm_real(flm_rec, f_cur, f_scal, &parameters);
    time_end = clock();
    printf("  - curvelet synthesis  : %4.4f seconds\n",
           (time_end - time_start) / (double)CLOCKS_PER_SEC);
    
    // Compute the maximum absolute error on the harmonic coefficients
    // printf("  -analysis, flm[1] :  %e+%ei\n",flm[1]);
    // printf("  -synthesis, flm_rec[1] :  %e+%ei\n",flm_rec[1]);
    printf("  - Maximum abs error  : %6.5e\n",
           maxerr_cplx(flm, flm_rec, L*L));fflush(NULL);
    
    free(flm);
    free(flm_rec);
    free(f_cur);
    free(f_scal);
}



/*!
 * Test that the directional algorithms reduce to the axisymmetric ones,
 * provided spin = 0, N = 1
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_transform_axisym_vs_directional_mw_test(B, L, J_min, seed)
{
    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = 1;
    parameters.upsample = 1;

    int spin = parameters.spin = 0;
    int J = s2let_j_max(&parameters);
    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;

    parameters.normalization = S2LET_WAV_NORM_DEFAULT;
    parameters.original_spin = 0;

    double wav_error, scal_error;

    complex double *f, *flm;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_mw(&f, L);

    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm(flm, L, spin, seed);

    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

    // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
    // from both transforms.
    complex double *f_wav_axisym, *f_scal_axisym, *f_wav_dir, *f_scal_dir;
    s2let_transform_axisym_allocate_mw_f_wav(&f_wav_axisym, &f_scal_axisym, &parameters);
    s2let_allocate_f_wav(&f_wav_dir, &f_scal_dir, &parameters);

    // Do both transforms
    s2let_transform_axisym_wav_analysis_mw(f_wav_axisym, f_scal_axisym, f, &parameters);
    s2let_analysis_px2wav(f_wav_dir, f_scal_dir, f, &parameters);

    // Compute the maximum absolute error in the computed wavelet transform
    wav_error = maxerr_cplx(f_wav_axisym, f_wav_dir, (J-J_min+1)*L*(2*L-1));
    scal_error = maxerr_cplx(f_scal_axisym, f_scal_dir, L*(2*L-1));

    printf("  - Maximum abs error in wavelets :         %6.5e\n", wav_error);
    printf("  - Maximum abs error in scaling function : %6.5e\n", scal_error);
    fflush(NULL);
}

/*!
 * Test that the directional multi-resolution algorithms reduce to
 * the axisymmetric ones, provided spin = 0, N = 1
 *
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  seed Random seed.
 * \retval none
 */
void s2let_transform_axisym_vs_directional_mw_multires_test(B, L, J_min, seed)
{
    s2let_parameters_t parameters = {};
    parameters.B = B;
    parameters.L = L;
    parameters.J_min = J_min;
    parameters.N = 1;

    int spin = parameters.spin = 0;

    int J = s2let_j_max(&parameters);
    int verbosity = parameters.verbosity = 0;
    ssht_dl_method_t dl_method = parameters.dl_method = SSHT_DL_RISBO;

    parameters.normalization = S2LET_WAV_NORM_DEFAULT;
    parameters.original_spin = 0;

    int samples, bandlimit, j;

    double wav_error, scal_error;

    complex double *f, *flm;
    s2let_allocate_lm(&flm, L);
    s2let_allocate_mw(&f, L);

    // Generate random harmonic coefficients for a complex signal
    s2let_lm_random_flm(flm, L, spin, seed);

    // Construct the corresponding signal on the sphere (MW sampling)
    ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

    // Allocate space for wavelet maps on the sphere (corresponding to the triplet B/L/J_min)
    // from both transforms.
    complex double *f_wav_axisym, *f_scal_axisym, *f_wav_dir, *f_scal_dir;
    s2let_transform_axisym_allocate_mw_f_wav_multires(&f_wav_axisym, &f_scal_axisym, &parameters);
    s2let_allocate_f_wav(&f_wav_dir, &f_scal_dir, &parameters);

    // Do both transforms
    s2let_transform_axisym_wav_analysis_mw_multires(f_wav_axisym, f_scal_axisym, f, &parameters);
    s2let_analysis_px2wav(f_wav_dir, f_scal_dir, f, &parameters);

    samples = 0;
    for (j = J_min; j <= J; ++j)
    {
        bandlimit = MIN(s2let_bandlimit(j, &parameters), L);
        samples += bandlimit * (2 * bandlimit - 1);
    }

    // Compute the maximum absolute error in the computed wavelet transform
    wav_error = maxerr_cplx(f_wav_axisym, f_wav_dir, samples);
    bandlimit = MIN(s2let_bandlimit(J_min-1, &parameters), L);
    scal_error = maxerr_cplx(f_scal_axisym, f_scal_dir, bandlimit*(2*bandlimit-1));

    printf("  - Maximum abs error in wavelets :         %6.5e\n", wav_error);
    printf("  - Maximum abs error in scaling function : %6.5e\n", scal_error);
    fflush(NULL);
}




void s2let_transform_performance_test(int B, int J_min, int NREPEAT, int NSCALE, int seed)
{
  s2let_parameters_t parameters = {};
  parameters.B = B;
  parameters.J_min = J_min;

  complex double *f, *flm, *flm_rec, *f_rec, *f_wav, *f_scal;
  clock_t time_start, time_end;
  int sc, repeat;
  double tottime_analysis = 0, tottime_synthesis = 0;
  double accuracy = 0.0;

  int L = 2;

  for (sc=0; sc<NSCALE; sc++) {

    L *= 2;

    parameters.L = L;

    s2let_allocate_lm(&flm, L);

    printf(" > L =  %i \n", L);
    for (repeat=0; repeat<NREPEAT; repeat++){

      printf("  -> Iteration : %i on %i \n",repeat+1,NREPEAT);

      s2let_lm_random_flm(flm, L, 0, seed);
      s2let_allocate_mw(&f, L);
      s2let_mw_alm2map(f, flm, L);
      s2let_transform_axisym_allocate_mw_f_wav(&f_wav, &f_scal, &parameters);

      time_start = clock();
      s2let_transform_axisym_wav_analysis_mw(f_wav, f_scal, f, &parameters);
      time_end = clock();
      tottime_synthesis += (time_end - time_start) / (double)CLOCKS_PER_SEC;
      //printf("  - Duration for S2LET synthesis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);

      s2let_allocate_mw(&f_rec, L);

      time_start = clock();
      s2let_transform_axisym_wav_synthesis_mw(f_rec, f_wav, f_scal, &parameters);
      time_end = clock();
      tottime_analysis += (time_end - time_start) / (double)CLOCKS_PER_SEC;

      s2let_allocate_lm(&flm_rec, L);
      s2let_mw_map2alm(flm_rec, f_rec, L);

      //printf("  - Duration for S2LET analysis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);
      accuracy += maxerr_cplx(flm, flm_rec, L*L);fflush(NULL);

      free(f);
      free(f_rec);
      free(flm_rec);
      free(f_wav);
      free(f_scal);

    }

    tottime_synthesis = tottime_synthesis / (double)NREPEAT;
    tottime_analysis = tottime_analysis / (double)NREPEAT;
    accuracy = accuracy / (double)NREPEAT;

    printf("  - Average duration for S2LET synthesis  : %5.5f seconds\n", tottime_synthesis);
    printf("  - Average duration for S2LET analysis   : %5.5f seconds\n", tottime_analysis);
    printf("  - Average max error on reconstruction  : %6.5e\n", accuracy);

    free(flm);

  }

}



void s2let_transform_performance_multires_test(int B, int J_min, int NREPEAT, int NSCALE, int seed)
{
  s2let_parameters_t parameters = {};
  parameters.B = B;
  parameters.J_min = J_min;

  complex double *f, *flm, *flm_rec, *f_rec, *f_wav, *f_scal;
  clock_t time_start, time_end;
  int sc, repeat;
  double tottime_analysis = 0, tottime_synthesis = 0;
  double accuracy = 0.0;

  int L = 2;

  for (sc=0; sc<NSCALE; sc++) {

    L *= 2;

    parameters.L = L;

    s2let_allocate_lm(&flm, L);

    printf(" > L =  %i \n", L);
    for (repeat=0; repeat<NREPEAT; repeat++){

      printf("  -> Iteration : %i on %i \n",repeat+1,NREPEAT);

      s2let_lm_random_flm(flm, L, 0, seed);
      s2let_allocate_mw(&f, L);
      s2let_mw_alm2map(f, flm, L);
      s2let_transform_axisym_allocate_mw_f_wav_multires(&f_wav, &f_scal, &parameters);

      time_start = clock();
      s2let_transform_axisym_wav_analysis_mw_multires(f_wav, f_scal, f, &parameters);
      time_end = clock();
      tottime_synthesis += (time_end - time_start) / (double)CLOCKS_PER_SEC;
      //printf("  - Duration for S2LET synthesis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);

      s2let_allocate_mw(&f_rec, L);

      time_start = clock();
      s2let_transform_axisym_wav_synthesis_mw_multires(f_rec, f_wav, f_scal, &parameters);
      time_end = clock();
      tottime_analysis += (time_end - time_start) / (double)CLOCKS_PER_SEC;

      s2let_allocate_lm(&flm_rec, L);
      s2let_mw_map2alm(flm_rec, f_rec, L);

      //printf("  - Duration for S2LET analysis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);
      accuracy += maxerr_cplx(flm, flm_rec, L*L);

      free(f);
      free(f_rec);
      free(flm_rec);
      free(f_wav);
      free(f_scal);

    }

    tottime_synthesis = tottime_synthesis / (double)NREPEAT;
    tottime_analysis = tottime_analysis / (double)NREPEAT;
    accuracy = accuracy / (double)NREPEAT;

    printf("  - Average duration for S2LET multires synthesis  : %5.5f seconds\n", tottime_synthesis);
    printf("  - Average duration for S2LET multires analysis   : %5.5f seconds\n", tottime_analysis);
    printf("  - Average max error on reconstruction  : %6.5e\n", accuracy);

    free(flm);

  }

}




void s2let_transform_lm_performance_test(int B, int J_min, int NREPEAT, int NSCALE, int seed)
{
  s2let_parameters_t parameters = {};
  parameters.B = B;
  parameters.J_min = J_min;

  complex double *flm, *flm_rec, *f_wav_lm, *f_scal_lm;
  clock_t time_start, time_end;
  int sc, repeat;
  double tottime_analysis = 0, tottime_synthesis = 0;
  double accuracy = 0.0;
  double *wav_lm, *scal_lm;

  int L = 4;

  for (sc=0; sc<NSCALE; sc++) {

    L *= 2;

    parameters.L = L;

    s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, &parameters);
    s2let_transform_axisym_lm_wav(wav_lm, scal_lm, &parameters);
    s2let_allocate_lm(&flm, L);

    printf(" > L =  %i \n", L);
    for (repeat=0; repeat<NREPEAT; repeat++){

      //printf("  -> Iteration : %i on %i \n",repeat+1,NREPEAT);

      s2let_lm_random_flm(flm, L, 0, seed);

      s2let_transform_axisym_lm_allocate_f_wav(&f_wav_lm, &f_scal_lm, &parameters);

      time_start = clock();
      s2let_transform_axisym_lm_wav_analysis(f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, &parameters);
      time_end = clock();
      tottime_synthesis += (time_end - time_start) / (double)CLOCKS_PER_SEC;
      //printf("  - Duration for S2LET synthesis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);

      s2let_allocate_lm(&flm_rec, L);

      time_start = clock();
      s2let_transform_axisym_lm_wav_synthesis(flm_rec, f_wav_lm, f_scal_lm, wav_lm, scal_lm, &parameters);
      time_end = clock();
      tottime_analysis += (time_end - time_start) / (double)CLOCKS_PER_SEC;

      //printf("  - Duration for S2LET analysis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);
      accuracy += maxerr_cplx(flm, flm_rec, L*L);

      free(flm_rec);
      free(f_wav_lm);
      free(f_scal_lm);

    }

    tottime_synthesis = tottime_synthesis / (double)NREPEAT;
    tottime_analysis = tottime_analysis / (double)NREPEAT;
    accuracy = accuracy / (double)NREPEAT;

    //printf("  - Average duration for S2LET synthesis  : %5.2e seconds\n", tottime_synthesis);
    //printf("  - Average duration for S2LET analysis   : %5.2e seconds\n", tottime_analysis);
    printf("  - Average duration for S2LET           : %5.2e seconds\n", (tottime_analysis+tottime_synthesis)/2);
    printf("  - Average max error on reconstruction  : %6.5e\n", accuracy);

    free(flm);
    free(wav_lm);
    free(scal_lm);

  }

}



void s2let_transform_lm_performance_multires_test(int B, int J_min, int NREPEAT, int NSCALE, int seed)
{
  s2let_parameters_t parameters = {};
  parameters.B = B;
  parameters.J_min = J_min;

  complex double *flm, *flm_rec, *f_wav_lm, *f_scal_lm;
  clock_t time_start, time_end;
  int sc, repeat;
  double tottime_analysis = 0, tottime_synthesis = 0;
  double accuracy = 0.0;
  double *wav_lm, *scal_lm;

  int L = 4;

  for (sc=0; sc<NSCALE; sc++) {

    L *= 2;

    parameters.L = L;

    s2let_transform_axisym_lm_allocate_wav(&wav_lm, &scal_lm, &parameters);
    s2let_transform_axisym_lm_wav(wav_lm, scal_lm, &parameters);
    s2let_allocate_lm(&flm, L);

    printf(" > L =  %i \n", L);
    for (repeat=0; repeat<NREPEAT; repeat++){

      //printf("  -> Iteration : %i on %i \n",repeat+1,NREPEAT);

      s2let_lm_random_flm(flm, L, 0, seed);

      s2let_transform_axisym_lm_allocate_f_wav_multires(&f_wav_lm, &f_scal_lm, &parameters);

      time_start = clock();
      s2let_transform_axisym_lm_wav_analysis_multires(f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, &parameters);
      time_end = clock();
      tottime_synthesis += (time_end - time_start) / (double)CLOCKS_PER_SEC;
      //printf("  - Duration for S2LET synthesis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);

      s2let_allocate_lm(&flm_rec, L);

      time_start = clock();
      s2let_transform_axisym_lm_wav_synthesis_multires(flm_rec, f_wav_lm, f_scal_lm, wav_lm, scal_lm, &parameters);
      time_end = clock();
      tottime_analysis += (time_end - time_start) / (double)CLOCKS_PER_SEC;

      //printf("  - Duration for S2LET analysis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);
      accuracy += maxerr_cplx(flm, flm_rec, L*L);

      free(flm_rec);
      free(f_wav_lm);
      free(f_scal_lm);

    }

    tottime_synthesis = tottime_synthesis / (double)NREPEAT;
    tottime_analysis = tottime_analysis / (double)NREPEAT;
    accuracy = accuracy / (double)NREPEAT;

    //printf("  - Average duration for S2LET multires synthesis  : %5.2e seconds\n", tottime_synthesis);
    //printf("  - Average duration for S2LET multires analysis   : %5.2e seconds\n", tottime_analysis);
    printf("  - Average duration for S2LET multires  : %5.2e seconds\n", (tottime_analysis+tottime_synthesis)/2);
    printf("  - Average max error on reconstruction  : %6.5e\n", accuracy);

    free(flm);
    free(wav_lm);
    free(scal_lm);

  }

}


int main(int argc, char *argv[])
{
  const int L = 16;
  const int N = L;
  const int B = 2;
  const int J_min = 0;
  const int spin = 0;

  s2let_parameters_t parameters = {};

  parameters.B = B;
  parameters.L = L;
  parameters.J_min = J_min;
  parameters.N = N;
  parameters.spin = spin;

  // This is too often zero, so we add 1 (zero will result in all random
  // numbers being the same).
  const int seed = (int)((double)clock()/(double)CLOCKS_PER_SEC) + 1;
  int l_min = s2let_L0(J_min, &parameters);

  printf("===========================================================================\n");
  printf("Testing S2LET facilities with the MW sampling\n");
  printf("===========================================================================\n");
  printf("PARAMETERS: ");
  printf("L = %i  N = %i  B = %i  l_wav_min = %i  spin = %i  seed = %i\n",
         L, N, B, l_min, spin, seed);
  //s2let_switch_wavtype(3);
  printf("---------------------------------------------------------------------------\n");
  printf("> Testing logfact binomial coefficient implementation...\n");
  // Don't use more than 62 as the argument.
  s2let_binomial_coefficient_test(62);
  printf("===========================================================================\n");
  printf("> Testing axisymmetric kernels...\n");
  s2let_tiling_axisym_test(B, L, J_min);
  printf("---------------------------------------------------------------------------\n");
  printf("> Testing directionality components...\n");
  s2let_tiling_direction_test(L, N);
  printf("---------------------------------------------------------------------------\n");
  printf("> Testing directional wavelets...\n");
  s2let_tiling_wavelet_test(B, L, J_min, N, spin);
  printf("===========================================================================\n");
  printf("> Testing axisymmetric wavelets in harmonic space...\n");
  s2let_transform_axisym_lm_wav_test(B, L, J_min, seed);
  printf("---------------------------------------------------------------------------\n");
  printf("> Testing axisymmetric multiresolution algorithm in harmonic space...\n");
  s2let_transform_axisym_lm_wav_multires_test(B, L, J_min, seed);
  printf("---------------------------------------------------------------------------\n");
  printf("> Testing directional wavelets in harmonic space...\n");
  s2let_wav_transform_harmonic_test(B, L, J_min, N, spin, seed);
  printf("---------------------------------------------------------------------------\n");
  printf("> Testing directional multiresolution algorithm in harmonic space...\n");
  s2let_wav_transform_harmonic_multires_test(B, L, J_min, N, spin, seed);
  printf("===========================================================================\n");
  printf("> Testing axisymmetric wavelets in pixel space...\n");
  s2let_transform_axisym_wav_test(B, L, J_min, seed);
  printf("---------------------------------------------------------------------------\n");
  printf("> Testing axisymmetric multiresolution algorithm in pixel space...\n");
  s2let_transform_axisym_wav_multires_test(B, L, J_min, seed);
  printf("---------------------------------------------------------------------------\n");
  printf("> Testing directional wavelets in pixel space...\n");
  s2let_wav_transform_mw_test(B, L, J_min, N, spin, seed);
  printf("---------------------------------------------------------------------------\n");
  printf("> Testing directional multiresolution algorithm in pixel space...\n");
  s2let_wav_transform_mw_multires_test(B, L, J_min, N, spin, seed);
  printf("===========================================================================\n");
  printf("> Comparing directional and axisymmetric algorithm in pixel space...\n");
  s2let_transform_axisym_vs_directional_mw_test(B, L, J_min, seed);
  printf("---------------------------------------------------------------------------\n");
  printf("> Comparing directional and axisymmetric multiresolution algorithm\n");
  printf("  in pixel space...\n");
  s2let_transform_axisym_vs_directional_mw_multires_test(B, L, J_min, seed);
  printf("===========================================================================\n");
  printf("> Testing real axisymmetric wavelets in pixel space...\n");
  s2let_transform_axisym_wav_real_test(B, L, J_min, seed);
  printf("---------------------------------------------------------------------------\n");
  printf("> Testing real axisymmetric multiresolution algorithm...\n");
  s2let_transform_axisym_wav_multires_real_test(B, L, J_min, seed);
  printf("---------------------------------------------------------------------------\n");
  printf("> Testing real directional wavelets in pixel space...\n");
  s2let_wav_transform_mw_real_test(B, L, J_min, N, seed);
  printf("---------------------------------------------------------------------------\n");
  printf("> Testing real directional multiresolution algorithm...\n");
  s2let_wav_transform_mw_multires_real_test(B, L, J_min, N, seed);
  printf("===========================================================================\n");
  printf("> Testing directional harmonic-to-wavelet transform...\n");
  s2let_wav_transform_lm2wav_test(B, L, J_min, N, spin, seed);
  printf("---------------------------------------------------------------------------\n");
  printf("> Testing directional multiresolution harmonic-to-wavelet transform...\n");
  s2let_wav_transform_lm2wav_multires_test(B, L, J_min, N, spin, seed);
  printf("---------------------------------------------------------------------------\n");
  printf("> Testing real directional harmonic-to-wavelet transform...\n");
  s2let_wav_transform_lm2wav_real_test(B, L, J_min, N, seed);
  printf("---------------------------------------------------------------------------\n");
  printf("> Testing real directional multiresolution harmonic-to-wavelet transform...\n");
  s2let_wav_transform_lm2wav_multires_real_test(B, L, J_min, N, seed);
  printf("===========================================================================\n");
  printf("> Testing directional wavelet transform with MWSS...\n");
  s2let_wav_transform_mwss_test(B, L, J_min, N, spin, seed);
  printf("---------------------------------------------------------------------------\n");
  printf("> Testing directional multiresolution wavelet transform with MWSS...\n");
  s2let_wav_transform_mwss_multires_test(B, L, J_min, N, spin, seed);
  printf("---------------------------------------------------------------------------\n");
  printf("> Testing real directional wavelet transform with MWSS...\n");
  s2let_wav_transform_mwss_real_test(B, L, J_min, N, seed);
  printf("---------------------------------------------------------------------------\n");
  printf("> Testing real directional multiresolution wavelet transform with MWSS...\n");
  s2let_wav_transform_mwss_multires_real_test(B, L, J_min, N, seed);
  printf("---------------------------------------------------------------------------\n");
  printf("---------------------------------------------------------------------------\n");
  printf("--------------------------- CURVELETS -------------------------------------\n");
  printf("===========================================================================\n");
  printf("---------------------------------------------------------------------------\n");
    printf("> Testing directionality of curvelets (i.e. slm) ...\n");
    s2let_tiling_curvelet_direction_test(L, N);
    printf("---------------------------------------------------------------------------\n");
    printf("> Testing curvelets...\n");
    s2let_tiling_curvelet_test(B, L, J_min, N, spin);
    printf("---------------------------------------------------------------------------\n");
    printf("> Testing curvelets in harmonic space...(s2let_cur_transform_harmonic_test)\n");
    s2let_cur_transform_harmonic_test(B, L, J_min, N, spin, seed);
    printf("---------------------------------------------------------------------------\n");
    printf("> Testing curvelet multiresolution algorithm in harmonic space...(s2let_cur_transform_harmonic_multires_test)\n");
    s2let_cur_transform_harmonic_multires_test(B, L, J_min, N, spin, seed);
    printf("===========================================================================\n");
    
    printf("---------------------------------------------------------------------------\n");
    printf("> Testing real curvelets in pixel space...\n");
    s2let_cur_transform_mw_real_test(B, L, J_min, N, seed);
    printf("---------------------------------------------------------------------------\n");
    printf("> Testing real multiresolution algorithm...\n");
    s2let_cur_transform_mw_multires_real_test(B, L, J_min, N, seed);
    printf("===========================================================================\n");
    printf("> Testing harmonic-to-curvelet transform...(s2let_cur_transform_lm2cur_test)\n");
    s2let_cur_transform_lm2cur_test(B, L, J_min, N, spin, seed);
    printf("---------------------------------------------------------------------------\n");
    printf("> Testing  multiresolution harmonic-to-curvelet transform...\n");
    s2let_cur_transform_lm2cur_multires_test(B, L, J_min, N, spin, seed);
    printf("---------------------------------------------------------------------------\n");
    printf("> Testing real harmonic-to-curvelet transform...\n");
    s2let_cur_transform_lm2cur_real_test(B, L, J_min, N, seed);
    printf("---------------------------------------------------------------------------\n");
    printf("> Testing real multiresolution harmonic-to-curvelet transform...\n");
    s2let_cur_transform_lm2cur_multires_real_test(B, L, J_min, N, seed);
    
    printf("===========================================================================\n");
    
    printf("---------------------------------------------------------------------------\n");
    printf("> Testing  curvelets in pixel space...(s2let_cur_transform_mw_test)\n");
    s2let_cur_transform_mw_test(B, L, J_min, N, spin, seed);
    printf("---------------------------------------------------------------------------\n");
    printf("> Testing curvelets multiresolution algorithm in pixel space...(s2let_cur_transform_mw_multires_test)\n");
    s2let_cur_transform_mw_multires_test(B, L, J_min, N, spin, seed);
    //  printf("===========================================================================\n");
    //  printf("> Comparing curvelet and axisymmetric algorithm in pixel space...\n");
    //  s2let_transform_axisym_vs_curvelet_mw_test(B, L, J_min, seed);
    
    // Do both transforms
    //  s2let_transform_axisym_cur_analysis_mw(f_cur_axisym, f_scal_axisym, f, &parameters);
    //  s2let_analysis_px2cur(f_cur_dir, f_scal_dir, f, &parameters);
    
    // Compute the maximum absolute error in the computed wavelet transform
    // cur_error = maxerr_cplx(f_cur_axisym, f_cur_dir, (J-J_min+1)*L*(2*L-1));
    // scal_error = maxerr_cplx(f_scal_axisym, f_scal_dir, L*(2*L-1));
    
    //  printf("---------------------------------------------------------------------------\n");
    //  printf("> Comparing curvelet and axisymmetric multiresolution algorithm\n");
    //  printf("  in pixel space...\n");
    //  s2let_transform_axisym_vs_curvelet_mw_multires_test(B, L, J_min, seed);

    printf("===========================================================================\n");
    printf("> Testing curvelet transform with MWSS...(s2let_cur_transform_mwss_test)\n");
    s2let_cur_transform_mwss_test(B, L, J_min, N, spin, seed);
    printf("---------------------------------------------------------------------------\n");
    printf("> Testing multiresolution curvelet transform with MWSS...\n");
    s2let_cur_transform_mwss_multires_test(B, L, J_min, N, spin, seed);
    printf("---------------------------------------------------------------------------\n");
    printf("> Testing real curvelet transform with MWSS...\n");
    s2let_cur_transform_mwss_real_test(B, L, J_min, N, seed);
    printf("---------------------------------------------------------------------------\n");
    printf("> Testing real multiresolution curvelet transform with MWSS...\n");
    s2let_cur_transform_mwss_multires_real_test(B, L, J_min, N, seed);


    /*
  const int NREPEAT = 50;
  const int NSCALE = 9;
    printf("==========================================================\n");
    printf("> Testing performances at full resolution...\n");
    s2let_transform_lm_performance_test(B, J_min, NREPEAT, NSCALE, seed);
    printf("----------------------------------------------------------\n");
    printf("> Testing performances with multiresolution...\n");
    s2let_transform_lm_performance_multires_test(B, J_min, NREPEAT, NSCALE, seed);
    printf("----------------------------------------------------------\n");
    printf("> Testing performances at full resolution...\n");
    //s2let_transform_performance_test(B, J_min, NREPEAT, NSCALE, seed);
    printf("----------------------------------------------------------\n");
    printf("> Testing performances with multiresolution...\n");
    //s2let_transform_performance_multires_test(B, J_min, NREPEAT, NSCALE, seed);
    */
  printf("===========================================================================\n");



  return 0;
}
