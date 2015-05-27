// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen
// CUREVELETS: added by Jennifer Y.H. Chan

#include "s2let.h"
#include <ssht.h>
#include <assert.h>
#include <complex.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*!
 * Test the identity relation of the directionality components for
 * the curvelet tiling in harmonic space.
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
    // ensuring exactness of the curvelet transform.
    error = s2let_tiling_direction_check_identity(s_elm, &parameters);
    printf("  - Maximum error : %6.5e\n", error);
    
    free(s_elm);
}

/*!
 * Test the identity relation of the curvelets for
 * the curvelet tiling in harmonic space.
 *
 * \param[in]  B Curvelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First curvelet scale to be used.
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
    // ensuring exactness of the curvelet transform.
    error = s2let_tiling_curvelet_check_identity(phi, psi, &parameters);
    printf("  - Maximum error : %6.5e\n", error);
    
    free(phi);
    free(psi);
}


/*!
 * Test the exactness of the full resolution curvelet transform
 * in harmonic space.
 *
 * \param[in]  B Curvelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First curvelet scale to be used.
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
 * \param[in]  B Curvelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First curvelet scale to be used.
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
 * Test the exactness of the full resolution curvelet transform
 * in pixel space for complex functions.
 *
 * \param[in]  B Curvelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First curvelet scale to be used.
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
 * \param[in]  B Curvelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First curvelet scale to be used.
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


/* ===========================s2let_test ==================================*/
// Define test parameters:
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
  printf("L = %i  N = %i  B = %i  l_cur_min = %i  spin = %i  seed = %i\n",
         L, N, B, l_min, spin, seed);
  printf("--------------------------- CURVELETS -------------------------------------\n");
  printf("===========================================================================\n");
  printf("---------------------------------------------------------------------------\n");
    printf("> Testing directionality of curvelets (i.e. slm) ...\n");
    s2let_tiling_curvelet_direction_test(L, N);
    printf("---------------------------------------------------------------------------\n");
    printf("> Testing curvelets...\n");
    s2let_tiling_curvelet_test(B, L, J_min, N, spin);
    printf("---------------------------------------------------------------------------\n");
    printf("> Testing curvelets in harmonic space...(s2let_cur_transform_harmonic_test: ana_cur_lm2lmn VS syn_cur_lmn2lm )\n");
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
    
    // Compute the maximum absolute error in the computed curvelet transform
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
  printf("===========================================================================\n");



  return 0;
}
