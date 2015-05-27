// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <so3.h>
#include <ssht.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/*!
 * Wavelet synthesis from Wigner space to harmonic space for complex signals.
 *
 * \param[out]  flm Spherical harmonic coefficients of input function.
 * \param[in]  f_wav_lmn Wavelet transform (Wigner coefficients of wavelet contribution).
 * \param[in]  f_scal_lm Wavelet transform (spherical harmonic coefficients of scaling contribution).
 * \param[in]  wav_lm Wavelet kernels.
 * \param[in]  scal_l Scaling function kernels.
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        s2let_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link s2let_synthesis_lmn2lm_real
 *                        \endlink instead for real signals.
 * \retval none
 */
void s2let_synthesis_lmn2lm(
    complex double *flm,
    const complex double *f_wav_lmn,
    const complex double *f_scal_lm,
    const complex double *wav_lm,
    const double *scal_l,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;
    int spin = parameters->spin;

    int j, el, m ,n;
    int J = s2let_j_max(parameters);
    int bandlimit = L;
    int Nj = N;

    int lm_ind, lmn_ind;
    so3_parameters_t so3_parameters = {};
    fill_so3_parameters(&so3_parameters, parameters);

    complex double psi;
    double phi;

    
    // For debugging:
    // Open data file '"f_wav_lmn.dat"' to write out f_wav_lm
    FILE *fp7, *fp8, *fp9, *fp10,*fp11, *fp12;
    fp7=fopen("4a_f_wav_lmn_syn.dat", "w");
    fp8=fopen("4a_f_wav_lmnONLYj0_syn.dat", "w");
    fp9=fopen("4a_f_wav_lmnONLYj1_syn.dat", "w");
    fp10=fopen("4a_f_wav_lmnONLYj2_syn.dat", "w");
    fp11=fopen("4a_f_wav_lmnONLYj3_syn.dat", "w");
    fp12=fopen("4a_f_wav_lmnONLYj4_syn.dat", "w");
    
    // For debugging:
    // Open data file '"f_lm.dat"' to write out f_lm
    FILE *fp, *fp2, *fp3, *fp4,*fp5, *fp6;
    fp=fopen("4b_f_lm_wav_syn.dat", "w");
    fp2=fopen("4b_f_lm_wavONLYj0_syn.dat", "w");
    fp3=fopen("4b_f_lm_wavONLYj1_syn.dat", "w");
    fp4=fopen("4b_f_lm_wavONLYj2_syn.dat", "w");
    fp5=fopen("4b_f_lm_wavONLYj3_syn.dat", "w");
    fp6=fopen("4b_f_lm_wavONLYj4_syn.dat", "w");
    
    
    int offset = 0;

    // Clear output
    for (el = 0; el < L*L; ++el)
        flm[el] = 0;

    for (j = J_min; j <= J; ++j)
    {
        if (!parameters->upsample)
        {
            bandlimit = MIN(s2let_bandlimit(j, parameters), L);
            so3_parameters.L = bandlimit;
            int Nj = MIN(N,bandlimit);
            Nj += (Nj+N)%2; // ensure N and Nj are both even or both odd
            so3_parameters.N = Nj;
        }

        for (n = -Nj+1; n < Nj; n+=1) //original n+=2
        {
            for (el = MAX(ABS(spin), ABS(n)); el < bandlimit; ++el)
            {
                ssht_sampling_elm2ind(&lm_ind, el, n);
                psi = wav_lm[j*L*L + lm_ind];
                for (m = -el; m <= el; ++m)
                {
                    ssht_sampling_elm2ind(&lm_ind, el, m);
                    so3_sampling_elmn2ind(&lmn_ind, el, m, n, &so3_parameters);
                    flm[lm_ind] += f_wav_lmn[offset + lmn_ind] * psi;
                    
                    // Write out to data file 'f_wav_lmnONLYj*_syn.dat'
                    fprintf(fp7, "%d, %d, %d, %d, %d, %d, %f, %f\n",j,n,el,m,offset,lmn_ind, creal(f_wav_lmn[offset +lmn_ind]), cimag(f_wav_lmn[offset +lmn_ind]));
                    if (j==J_min)
                        fprintf(fp8, "%f, %f\n", creal(f_wav_lmn[offset +lmn_ind]), cimag(f_wav_lmn[offset +lmn_ind]));
                    if (j==1)
                        fprintf(fp9, "%f, %f\n", creal(f_wav_lmn[offset +lmn_ind]), cimag(f_wav_lmn[offset +lmn_ind]));
                    if (j==2)
                        fprintf(fp10, "%f, %f\n", creal(f_wav_lmn[offset +lmn_ind]), cimag(f_wav_lmn[offset +lmn_ind]));
                    if (j==3)
                        fprintf(fp11, "%f, %f\n", creal(f_wav_lmn[offset +lmn_ind]), cimag(f_wav_lmn[offset +lmn_ind]));
                    if (j==4)
                        fprintf(fp12, "%f, %f\n", creal(f_wav_lmn[offset +lmn_ind]), cimag(f_wav_lmn[offset +lmn_ind]));
                    // Write out to data file 'f_lm_ONLYj*.dat_syn.dat'
                    fprintf(fp, "%d, %d, %d, %d, %d, %f, %f\n",j,n,el,m,lm_ind,  creal(flm[lm_ind]), cimag(flm[lm_ind]));
                    if (j==J_min)
                        fprintf(fp2, "%f, %f\n", creal(flm[lm_ind]), cimag(flm[lm_ind]));
                    if (j==1)
                        fprintf(fp3, "%f, %f\n", creal(flm[lm_ind]), cimag(flm[lm_ind]));
                    if (j==2)
                        fprintf(fp4, "%f, %f\n", creal(flm[lm_ind]), cimag(flm[lm_ind]));
                    if (j==3)
                        fprintf(fp5, "%f, %f\n", creal(flm[lm_ind]), cimag(flm[lm_ind]));
                    if (j==4)
                        fprintf(fp6, "%f, %f\n", creal(flm[lm_ind]), cimag(flm[lm_ind]));
                    
                }
            }
        }
        offset += so3_sampling_flmn_size(&so3_parameters);
    }
    
    if (!parameters->upsample)
        bandlimit = MIN(s2let_bandlimit(J_min-1, parameters), L);

    for (el = ABS(spin); el < bandlimit; ++el)
    {
        phi = sqrt(4.0*PI/(2*el+1)) * scal_l[el];
        for (m = -el; m <= el; ++m)
        {
            ssht_sampling_elm2ind(&lm_ind, el, m);
            flm[lm_ind] += f_scal_lm[lm_ind] * phi;
        }
    }
    
    
    // Close file '"flm.dat"'
    fclose(fp);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);
    fclose(fp5);
    fclose(fp6);
    fclose(fp7);
    fclose(fp8);
    fclose(fp9);
    fclose(fp10);
    fclose(fp11);
    fclose(fp12);

}

/*!
 * Wavelet synthesis from Wigner space to harmonic space for real signals.
 *
 * \param[out]  flm Spherical harmonic coefficients of input function.
 * \param[in]  f_wav_lmn Wavelet transform (Wigner coefficients of wavelet contribution).
 * \param[in]  f_scal_lm Wavelet transform (spherical harmonic coefficients of scaling contribution).
 * \param[in]  wav_lm Wavelet kernels.
 * \param[in]  scal_l Scaling function kernels.
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        s2let_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link s2let_synthesis_lmn2lm
 *                        \endlink instead for complex signals.
 * \retval none
 */
void s2let_synthesis_lmn2lm_real(
    complex double *flm,
    const complex double *f_wav_lmn,
    const complex double *f_scal_lm,
    const complex double *wav_lm,
    const double *scal_l,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;

    int j, el, m ,n;
    int J = s2let_j_max(parameters);
    int bandlimit = L;
    int Nj = N;

    int lm_ind, lmn_ind;
    so3_parameters_t so3_parameters = {};
    fill_so3_parameters(&so3_parameters, parameters);

    complex double psi, npsi;
    double phi;

    int offset = 0;

    // Clear output
    for (el = 0; el < L*L; ++el)
        flm[el] = 0;

    for (j = J_min; j <= J; ++j)
    {
        if (!parameters->upsample)
        {
            bandlimit = MIN(s2let_bandlimit(j, parameters), L);
            so3_parameters.L = bandlimit;
            int Nj = MIN(N,bandlimit);
            Nj += (Nj+N)%2; // ensure N and Nj are both even or both odd
            so3_parameters.N = Nj;
        }

        for (n = 1-Nj%2; n < Nj; n+=1)   //original n+=2
        {
            for (el = n; el < bandlimit; ++el)
            {
                ssht_sampling_elm2ind(&lm_ind, el, n);
                psi = wav_lm[j*L*L + lm_ind];

                if (n)
                {
                    ssht_sampling_elm2ind(&lm_ind, el, -n);
                    npsi = wav_lm[j*L*L + lm_ind];
                }

                for (m = -el; m <= el; ++m)
                {
                    ssht_sampling_elm2ind(&lm_ind, el, m);
                    so3_sampling_elmn2ind_real(&lmn_ind, el, m, n, &so3_parameters);
                    flm[lm_ind] += f_wav_lmn[offset + lmn_ind] * psi;

                    if (n)
                    {
                        so3_sampling_elmn2ind_real(&lmn_ind, el, -m, n, &so3_parameters);
                        int sign = (m+n)%2 ? -1 : 1;
                        flm[lm_ind] += sign * conj(f_wav_lmn[offset + lmn_ind]) * npsi;
                    }
                }
            }
        }
        offset += so3_sampling_flmn_size(&so3_parameters);
    }

    if (!parameters->upsample)
        bandlimit = MIN(s2let_bandlimit(J_min-1, parameters), L);

    for (el = 0; el < bandlimit; ++el)
    {
        phi = sqrt(4.0*PI/(2*el+1)) * scal_l[el];
        for (m = -el; m <= el; ++m)
        {
            ssht_sampling_elm2ind(&lm_ind, el, m);
            flm[lm_ind] += f_scal_lm[lm_ind] * phi;
        }
    }
}

/*!
 * Wavelet synthesis from wavelet space to harmonic space for complex signals.
 *
 * \param[out]  flm Spherical harmonic coefficients of the signal
 * \param[in]  f_wav Array of wavelets maps
 * \param[in]  f_scal Scaling function map
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        s2let_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link s2let_synthesis_wav2lm_real
 *                        \endlink instead for real signals.
 * \retval none
 */
void s2let_synthesis_wav2lm(
    complex double *flm,
    const complex double *f_wav,
    const complex double *f_scal,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;
    ssht_dl_method_t dl_method = parameters->dl_method;

    int bandlimit = L;
    int verbosity = 0;
    so3_parameters_t so3_parameters = {};
    fill_so3_parameters(&so3_parameters, parameters);

    int j, offset, offset_lmn;
    int J = s2let_j_max(parameters);

    complex double *wav_lm;
    double *scal_l;
    s2let_tiling_wavelet_allocate(&wav_lm, &scal_l, parameters);
    s2let_tiling_wavelet(wav_lm, scal_l, parameters);

    complex double *f_wav_lmn, *f_scal_lm;
    s2let_allocate_lmn_f_wav(&f_wav_lmn, &f_scal_lm, parameters);    

    if (!parameters->upsample)
        bandlimit = MIN(s2let_bandlimit(J_min-1, parameters), L);

    // Note, this is a spin-0 transform!
    switch (parameters->sampling_scheme)
    {
    case S2LET_SAMPLING_MW:
        ssht_core_mw_forward_sov_conv_sym(f_scal_lm, f_scal, bandlimit, 0, dl_method, verbosity);
        break;
    case S2LET_SAMPLING_MW_SS:
        ssht_core_mw_forward_sov_conv_sym_ss(f_scal_lm, f_scal, bandlimit, 0, dl_method, verbosity);
        break;
    default:
        S2LET_ERROR_GENERIC("Sampling scheme not supported.");
    }
    
    //printf(" ***** - size of f_wav :  %ld\n",sizeof(f_wav));
    //printf(" ***** - size of f_wav_lmn :  %ld\n",sizeof(f_wav_lmn));
    
    // For debugging:
    int arrayind_min,arrayind, arrayind_max;
    FILE *fp10;
    arrayind_min = 0;
    arrayind_max = (2*L-1)*L*(2*N-1);    //L*L;
    fp10=fopen("4aaa_check_f_wav.dat", "w");
    for (arrayind = arrayind_min; arrayind <  arrayind_max; arrayind++ )
    {
        fprintf(fp10, "%d, %f, %f\n", arrayind, creal(f_wav[arrayind]), cimag(f_wav[arrayind]));
    }
    fclose(fp10);
    
    
    
    // For debugging:
    // Open data file '"f_wav_lmn.dat"' to write out f_wav_lm    
    FILE *fp ;  //, *fp2, *fp3, *fp4,*fp5, *fp6;
    FILE *fp7; //*fp8, *fp9, *fp10,*fp11, *fp12, *fp13;
    FILE *fp13 ; //, *fp14;
    fp=fopen("4ab_f_wav_lmn_syn_wav2lm_so3coreforwardviassht.dat", "w");
    fp7=fopen("4aa_f_wav_syn_wav2lm_so3coreforwardviassht.dat", "w");
    fp13=fopen("4abc_wav_so3para.dat", "w");
    
    
    
    offset = 0;
    offset_lmn = 0;
    for (j = J_min; j <= J; ++j)
    {
        if (!parameters->upsample)
        {
            bandlimit = MIN(s2let_bandlimit(j, parameters), L);
            so3_parameters.L = bandlimit;
            int Nj = MIN(N,bandlimit);
            Nj += (Nj+N)%2; // ensure N and Nj are both even or both odd
            so3_parameters.N = Nj;
        }

        so3_parameters.L0 = s2let_L0(j, parameters);
        
        
        // For debugging
        // Write out to data file '"4abc_wav_so3para.dat"'
        fprintf(fp13, "%d,%d,%d,%d\n",j, so3_parameters.L0, so3_parameters.L, so3_parameters.N);
        // Write out to data file '"4ab_f_wav_lmn_syn_wav2lm_so3coreforwardviassht.dat"'
        //fprintf(fp, "%d, %f, %f, %ld, %d, %d, %f, %f, %ld\n",j ,creal(*f_wav_lmn), cimag(*f_wav_lmn) , sizeof(f_wav_lmn), offset_lmn, so3_sampling_flmn_size(&so3_parameters), creal(*f_wav_lmn +offset_lmn),cimag(*f_wav_lmn +offset_lmn), sizeof(f_wav_lmn +offset_lmn));
        fprintf(fp, "%d,  %6.5e+i%6.5e,  %ld, %d, %d,  %6.5e+i%6.5e,  %ld\n",j ,*(f_wav_lmn) , sizeof(f_wav_lmn), offset_lmn, so3_sampling_flmn_size(&so3_parameters), (*f_wav_lmn +offset_lmn), sizeof(f_wav_lmn +offset_lmn));
        // Write out to data file '"4aa_f_wav_syn_wav2lm_so3coreforwardviassht.dat"'
        fprintf(fp7, "%d,  %6.5e+i%6.5e,  %ld, %d,  %6.5e+i%6.5e, , %ld\n",j, *(f_wav), sizeof(f_wav), offset, (*f_wav+offset), sizeof(f_wav+offset));


        so3_core_forward_via_ssht(
            f_wav_lmn + offset_lmn,
            f_wav + offset,
            &so3_parameters
        );

        offset_lmn += so3_sampling_flmn_size(&so3_parameters);
        offset += so3_sampling_f_size(&so3_parameters);
        
    }

    // Close file:
    fclose(fp);
    fclose(fp7);
    fclose(fp13);
    
    s2let_synthesis_lmn2lm(flm, f_wav_lmn, f_scal_lm, wav_lm, scal_l, parameters);
    
    // For debugging:
    int arrayindex;
    FILE *fp9;
    fp9=fopen("4ac_check_f_wav_lmn.dat", "w");
    for (arrayindex =0; arrayindex < arrayind_max; arrayindex++ )
    {
        fprintf(fp9, "%d, %f, %f\n", arrayindex, creal(f_wav_lmn[arrayindex]), cimag(f_wav_lmn[arrayindex]));
    }
    fclose(fp9);
    
    // For debugging:
    arrayindex=0;
    FILE *fp11;
    fp11=fopen("4ad_check_flm_out_fromwav2lm.dat", "w");
    for (arrayindex =0; arrayindex < arrayind_max; arrayindex++ )
    {
        fprintf(fp11, "%d, %f, %f\n", arrayindex, creal(flm[arrayindex]), cimag(flm[arrayindex]));
    }
    fclose(fp11);

    free(wav_lm);
    free(scal_l);
    free(f_scal_lm);
    free(f_wav_lmn);

    
}

/*!
 * Wavelet synthesis from wavelet space to harmonic space for real signals.
 *
 * \param[out]  flm Spherical harmonic coefficients of the signal
 * \param[in]  f_wav Array of wavelets maps
 * \param[in]  f_scal Scaling function map
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        s2let_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link s2let_synthesis_wav2lm
 *                        \endlink instead for complex signals.
 * \retval none
 */
void s2let_synthesis_wav2lm_real(
    complex double *flm,
    const double *f_wav,
    const double *f_scal,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;
    ssht_dl_method_t dl_method = parameters->dl_method;

    s2let_parameters_t real_parameters = *parameters;
    real_parameters.reality = 1;

    int bandlimit = L;
    int verbosity = 0;
    so3_parameters_t so3_parameters = {};
    fill_so3_parameters(&so3_parameters, &real_parameters);

    int j, offset, offset_lmn;
    int J = s2let_j_max(&real_parameters);

    complex double *wav_lm;
    double *scal_l;
    s2let_tiling_wavelet_allocate(&wav_lm, &scal_l, &real_parameters);
    s2let_tiling_wavelet(wav_lm, scal_l, &real_parameters);

    complex double *f_wav_lmn, *f_scal_lm;
    s2let_allocate_lmn_f_wav(&f_wav_lmn, &f_scal_lm, &real_parameters);

    if (!parameters->upsample)
        bandlimit = MIN(s2let_bandlimit(J_min-1, &real_parameters), L);

    switch (parameters->sampling_scheme)
    {
    case S2LET_SAMPLING_MW:
        ssht_core_mw_forward_sov_conv_sym_real(f_scal_lm, f_scal, bandlimit, dl_method, verbosity);
        break;
    case S2LET_SAMPLING_MW_SS:
        ssht_core_mw_forward_sov_conv_sym_ss_real(f_scal_lm, f_scal, bandlimit, dl_method, verbosity);
        break;
    default:
        S2LET_ERROR_GENERIC("Sampling scheme not supported.");
    }


    offset = 0;
    offset_lmn = 0;
    for (j = J_min; j <= J; ++j)
    {
        if (!parameters->upsample)
        {
            bandlimit = MIN(s2let_bandlimit(j, &real_parameters), L);
            so3_parameters.L = bandlimit;
            int Nj = MIN(N,bandlimit);
            Nj += (Nj+N)%2; // ensure N and Nj are both even or both odd
            so3_parameters.N = Nj;
        }

        so3_parameters.L0 = s2let_L0(j, parameters);

        so3_core_forward_via_ssht_real(
            f_wav_lmn + offset_lmn,
            f_wav + offset,
            &so3_parameters
        );

        offset_lmn += so3_sampling_flmn_size(&so3_parameters);
        offset += so3_sampling_f_size(&so3_parameters);
    }

    s2let_synthesis_lmn2lm_real(flm, f_wav_lmn, f_scal_lm, wav_lm, scal_l, &real_parameters);

    free(wav_lm);
    free(scal_l);
    free(f_scal_lm);
    free(f_wav_lmn);
}

/*!
 * Wavelet synthesis from wavelet space to pixel space for complex signals.
 *
 * \param[out]  f Signal on the sphere
 * \param[in]  f_wav Array of wavelets maps
 * \param[in]  f_scal Scaling function map
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        s2let_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link s2let_synthesis_wav2px_real
 *                        \endlink instead for real signals.
 * \retval none
 */
void s2let_synthesis_wav2px(
    complex double *f,
    const complex double *f_wav,
    const complex double *f_scal,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int spin = parameters->spin;
    ssht_dl_method_t dl_method = parameters->dl_method;
    int verbosity = 0;

    complex double *flm;
    s2let_allocate_lm(&flm, L);

    s2let_synthesis_wav2lm(flm, f_wav, f_scal, parameters);

    switch (parameters->sampling_scheme)
    {
    case S2LET_SAMPLING_MW:
        ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);
        break;
    case S2LET_SAMPLING_MW_SS:
        ssht_core_mw_inverse_sov_sym_ss(f, flm, L, spin, dl_method, verbosity);
        break;
    default:
        S2LET_ERROR_GENERIC("Sampling scheme not supported.");
    }

    free(flm);
}

/*!
 * Wavelet synthesis from wavelet space to pixel space for real signals.
 *
 * \param[out]  f Signal on the sphere
 * \param[in]  f_wav Array of wavelets maps
 * \param[in]  f_scal Scaling function map
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        s2let_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link s2let_synthesis_wav2px
 *                        \endlink instead for complex signals.
 * \retval none
 */
void s2let_synthesis_wav2px_real(
    double *f,
    const double *f_wav,
    const double *f_scal,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    ssht_dl_method_t dl_method = parameters->dl_method;
    int verbosity = 0;

    complex double *flm;
    s2let_allocate_lm(&flm, L);

    s2let_synthesis_wav2lm_real(flm, f_wav, f_scal, parameters);

    switch (parameters->sampling_scheme)
    {
    case S2LET_SAMPLING_MW:
        ssht_core_mw_inverse_sov_sym_real(f, flm, L, dl_method, verbosity);
        break;
    case S2LET_SAMPLING_MW_SS:
        ssht_core_mw_inverse_sov_sym_ss_real(f, flm, L, dl_method, verbosity);
        break;
    default:
        S2LET_ERROR_GENERIC("Sampling scheme not supported.");
    }

    free(flm);
}
