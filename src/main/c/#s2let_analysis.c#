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
 * Wavelet analysis from harmonic space to Wigner space for complex signals.
 *
 * \param[out]  f_wav_lmn Wavelet transform (Wigner coefficients of wavelet contribution).
 * \param[out]  f_scal_lm Wavelet transform (Spherical harmonic coefficients of scaling contribution).
 * \param[in]  flm Spherical harmonic coefficients of input function.
 * \param[in]  wav_lm Wavelet kernels in harmonic space.
 * \param[in]  scal_l Scaling function kernels in harmonic space.
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        s2let_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link s2let_analysis_lm2lmn_real
 *                        \endlink instead for real signals.
 * \retval none
 */
void s2let_analysis_lm2lmn(
    complex double *f_wav_lmn,
    complex double *f_scal_lm,
    const complex double *flm,
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

    int offset = 0;

    for (j = J_min; j <= J; ++j)
    {

        if (!parameters->upsample)
        {
            bandlimit = MIN(s2let_bandlimit(j, parameters), L);
            so3_parameters.L = bandlimit;
            Nj = MIN(N,bandlimit);
            // ensure N and Nj are both even or both odd
            Nj += (Nj+N)%2;
            so3_parameters.N = Nj;
        }

        for (n = -Nj+1; n < Nj; n+=2)
        {
            for (el = MAX(ABS(spin), ABS(n)); el < bandlimit; ++el)
            {
                ssht_sampling_elm2ind(&lm_ind, el, n);
                psi = 8*PI*PI/(2*el+1) * conj(wav_lm[j*L*L + lm_ind]);
                for (m = -el; m <= el; ++m)
                {
                    ssht_sampling_elm2ind(&lm_ind, el, m);
                    so3_sampling_elmn2ind(&lmn_ind, el, m, n, &so3_parameters);
                    f_wav_lmn[offset + lmn_ind] = flm[lm_ind] * psi;
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
            f_scal_lm[lm_ind] = flm[lm_ind] * phi;
        }
    }
}

/*!
 * Wavelet analysis from harmonic space to Wigner space for real signals.
 *
 * \param[out]  f_wav_lmn Wavelet transform (Wigner coefficients of wavelet contribution).
 * \param[out]  f_scal_lm Wavelet transform (spherical harmonic coefficients of scaling contribution).
 * \param[in]  flm Spherical harmonic coefficients of input function.
 * \param[in]  wav_lm Wavelet kernels.
 * \param[in]  scal_l Scaling function kernels.
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        s2let_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link s2let_analysis_lm2lmn
 *                        \endlink instead for complex signals.
 * \retval none
 */
void s2let_analysis_lm2lmn_real(
    complex double *f_wav_lmn,
    complex double *f_scal_lm,
    const complex double *flm,
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

    complex double psi;
    double phi;

    int offset = 0;

    for (j = J_min; j <= J; ++j)
    {
        if (!parameters->upsample)
        {
            bandlimit = MIN(s2let_bandlimit(j, parameters), L);
            so3_parameters.L = bandlimit;
            int Nj = MIN(N,bandlimit);
            // ensure N and Nj are both even or both odd
            Nj += (Nj+N)%2;
            so3_parameters.N = Nj;
        }

        for (n = 1-Nj%2; n < Nj; n+=2)
        {
            for (el = n; el < bandlimit; ++el)
            {
                ssht_sampling_elm2ind(&lm_ind, el, n);
                psi = 8*PI*PI/(2*el+1) * conj(wav_lm[j*L*L + lm_ind]);
                for (m = -el; m <= el; ++m)
                {
                    ssht_sampling_elm2ind(&lm_ind, el, m);
                    so3_sampling_elmn2ind_real(&lmn_ind, el, m, n, &so3_parameters);
                    f_wav_lmn[offset + lmn_ind] = flm[lm_ind] * psi;
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
            f_scal_lm[lm_ind] = flm[lm_ind] * phi;
        }
    }
}


/*!
 * Wavelet analysis from harmonic space to wavelet space for complex signals.
 * with fully manual wavelet tiling, using multiresolution as default with 
 * the band-limits provided in input.
 *
 * \param[out]  f_wav Array of wavelet maps
 * \param[out]  f_scal Scaling function map
 * \param[in]  flm Spherical harmonic coefficients of the signal
 * \param[in]  scal_l Array of size L containing the \ell space tiling for the scaling fct.
                      It is only \ell because it is assumed to be axisymmetric.
 * \param[in]  wav_lm Array of size (J+1)*L*L containing the (\ell, space) harmonic coefs
                     of the wavelets. They can be directional. These must make sense and 
                     define a valid invertible transform as no extra checks are performed.
 * \param[in]  scal_bandlimit Same as wav_bandlimits but only one integer:
                      the band-limit of the scaling function.
 * \param[in]  wav_bandlimits Array of integers of size J+1 containing the band-limits 
                     of the wavelets. Will be used to do the multiresolution.
                     These must make sense and define a valid invertible transform
                     as no extra checks are performed.
 * \param[in]  J Number of scales in total (in wav_bandlimits) is J+1. 
 * \param[in]  L Band-limit for the transform: defines the size of all awways.
 * \param[in]  spin Spin (integer) to perform the transform 
 * \param[in]  N Azimuthal band-limit for the directional transform
 * \retval none
 */
void s2let_analysis_lm2wav_manual(
    complex double *f_wav,
    complex double *f_scal,
    const complex double *flm,
    const double *scal_l,
    const complex double *wav_lm,
    const int scal_bandlimit,
    const int *wav_bandlimits,
    int J,
    int L,
    int spin,
    int N
) {
    s2let_parameters_t parameters = {};
    parameters.L = L;
    parameters.J_min = 0;
    parameters.B = pow(L, 1.0/(float)(J+2));
    parameters.N = N;
    parameters.dl_method = SSHT_DL_RISBO;
    printf("J = %i\n", J);
    printf("B = %g\n", parameters.B);
    printf("N = %i\n", parameters.N);
    printf("L = %i\n", parameters.L);

    int bandlimit = L;
    int verbosity = 0;
    so3_parameters_t so3_parameters = {};
    fill_so3_parameters(&so3_parameters, &parameters);

    int j, offset, offset_lmn;
    complex double *f_wav_lmn, *f_scal_lm;
    complex double psi, npsi;
    double phi;
    int el, m, n, lm_ind, lmn_ind;

    int Nj = N;

    bandlimit = MIN(scal_bandlimit, L);
    f_scal_lm = (complex double*)calloc(bandlimit*bandlimit, sizeof(complex double));

    for (el = ABS(spin); el < bandlimit; ++el)
    {
        phi = sqrt(4.0*PI/(2*el+1)) * scal_l[el];
        for (m = -el; m <= el; ++m)
        {
            ssht_sampling_elm2ind(&lm_ind, el, m);
            f_scal_lm[lm_ind] = flm[lm_ind] * phi;
        }
    }

    switch (parameters.sampling_scheme)
    {
    case S2LET_SAMPLING_MW:
        ssht_core_mw_inverse_sov_sym(f_scal, f_scal_lm, bandlimit, 0, parameters.dl_method, verbosity);
        break;
    case S2LET_SAMPLING_MW_SS:
        ssht_core_mw_inverse_sov_sym_ss(f_scal, f_scal_lm, bandlimit, 0, parameters.dl_method, verbosity);
        break;
    default:
        S2LET_ERROR_GENERIC("Sampling scheme not supported.");
    }

    free(f_scal_lm);

    offset = 0;
    for (j = 0; j <= J; ++j)
    {
        bandlimit = MIN(wav_bandlimits[j], L);
        so3_parameters.L = bandlimit;
        int Nj = MIN(N, bandlimit);
        Nj += (Nj+N)%2; // ensure N and Nj are both even or both odd
        so3_parameters.N = Nj;
        so3_parameters.L0 = 0;

        f_wav_lmn = (complex double*)calloc(so3_sampling_flmn_size(&so3_parameters), sizeof(complex double));

        for (n = -Nj+1; n < Nj; n+=2)
        {
            for (el = MAX(ABS(spin), ABS(n)); el < bandlimit; ++el)
            {
                ssht_sampling_elm2ind(&lm_ind, el, n);
                psi = 8*PI*PI/(2*el+1) * conj(wav_lm[j*L*L + lm_ind]);
                for (m = -el; m <= el; ++m)
                {
                    ssht_sampling_elm2ind(&lm_ind, el, m);
                    so3_sampling_elmn2ind(&lmn_ind, el, m, n, &so3_parameters);
                    f_wav_lmn[lmn_ind] = flm[lm_ind] * psi;
                }
            }
        }

        so3_core_inverse_via_ssht(
            f_wav + offset,
            f_wav_lmn,
            &so3_parameters
        );

        free(f_wav_lmn);
        offset += so3_sampling_f_size(&so3_parameters);

    }

}


/*!
 * Wavelet analysis from harmonic space to wavelet space for complex signals.
 *
 * \param[out]  f_wav Array of wavelet maps
 * \param[out]  f_scal Scaling function map
 * \param[in]  flm Spherical harmonic coefficients of the signal
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        s2let_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link s2let_analysis_lm2wav_real
 *                        \endlink instead for real signals.
 * \retval none
 */
void s2let_analysis_lm2wav(
    complex double *f_wav,
    complex double *f_scal,
    const complex double *flm,
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
    s2let_analysis_lm2lmn(f_wav_lmn, f_scal_lm, flm, wav_lm, scal_l, parameters);

    if (!parameters->upsample)
        bandlimit = MIN(s2let_bandlimit(J_min-1, parameters), L);

    // Note, this is a spin-0 transform!
    switch (parameters->sampling_scheme)
    {
    case S2LET_SAMPLING_MW:
        ssht_core_mw_inverse_sov_sym(f_scal, f_scal_lm, bandlimit, 0, dl_method, verbosity);
        break;
    case S2LET_SAMPLING_MW_SS:
        ssht_core_mw_inverse_sov_sym_ss(f_scal, f_scal_lm, bandlimit, 0, dl_method, verbosity);
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
            bandlimit = MIN(s2let_bandlimit(j, parameters), L);
            so3_parameters.L = bandlimit;
            int Nj = MIN(N,bandlimit);
            Nj += (Nj+N)%2; // ensure N and Nj are both even or both odd
            so3_parameters.N = Nj;
        }

        so3_parameters.L0 = s2let_L0(j, parameters);

        so3_core_inverse_via_ssht(
            f_wav + offset,
            f_wav_lmn + offset_lmn,
            &so3_parameters
        );
        offset_lmn += so3_sampling_flmn_size(&so3_parameters);
        offset += so3_sampling_f_size(&so3_parameters);
    }

    free(wav_lm);
    free(scal_l);
    free(f_scal_lm);
    free(f_wav_lmn);
}

/*!
 * Wavelet analysis from harmonic space to wavlet space for real signals.
 *
 * \param[out]  f_wav Array of wavelet maps
 * \param[out]  f_scal Scaling function map
 * \param[in]  flm Spherical harmonic coefficients of the signal
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        s2let_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link s2let_analysis_lm2wav
 *                        \endlink instead for complex signals.
 * \retval none
 */
void s2let_analysis_lm2wav_real(
    double *f_wav,
    double *f_scal,
    const complex double *flm,
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
    s2let_analysis_lm2lmn_real(f_wav_lmn, f_scal_lm, flm, wav_lm, scal_l, &real_parameters);

    if (!parameters->upsample)
        bandlimit = MIN(s2let_bandlimit(J_min-1, &real_parameters), L);

    switch (parameters->sampling_scheme)
    {
    case S2LET_SAMPLING_MW:
        ssht_core_mw_inverse_sov_sym_real(f_scal, f_scal_lm, bandlimit, dl_method, verbosity);
        break;
    case S2LET_SAMPLING_MW_SS:
        ssht_core_mw_inverse_sov_sym_ss_real(f_scal, f_scal_lm, bandlimit, dl_method, verbosity);
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

        so3_core_inverse_via_ssht_real(
            f_wav + offset,
            f_wav_lmn + offset_lmn,
            &so3_parameters
        );
        offset_lmn += so3_sampling_flmn_size(&so3_parameters);
        offset += so3_sampling_f_size(&so3_parameters);
    }

    free(wav_lm);
    free(scal_l);
    free(f_scal_lm);
    free(f_wav_lmn);
}

/*!
 * Wavelet analysis from pixel space to wavelet space for complex signals.
 *
 * \param[out]  f_wav Array of wavelet maps
 * \param[out]  f_scal Scaling function map
 * \param[in]  f Signal on the sphere
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        s2let_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link s2let_analysis_px2wav_real
 *                        \endlink instead for real signals.
 * \retval none
 */
void s2let_analysis_px2wav(
    complex double *f_wav,
    complex double *f_scal,
    const complex double *f,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int spin = parameters->spin;
    ssht_dl_method_t dl_method = parameters->dl_method;
    int verbosity = parameters->verbosity;

    complex double *flm;
    s2let_allocate_lm(&flm, L);

    switch (parameters->sampling_scheme)
    {
    case S2LET_SAMPLING_MW:
        ssht_core_mw_forward_sov_conv_sym(flm, f, L, spin, dl_method, verbosity);
        break;
    case S2LET_SAMPLING_MW_SS:
        ssht_core_mw_forward_sov_conv_sym_ss(flm, f, L, spin, dl_method, verbosity);
        break;
    default:
        S2LET_ERROR_GENERIC("Sampling scheme not supported.");
    }

    s2let_analysis_lm2wav(f_wav, f_scal, flm, parameters);

    free(flm);
}

/*!
 * Wavelet analysis from pixel space to wavelet space for real signals.
 *
 * \param[out]  f_wav Array of wavelet maps
 * \param[out]  f_scal Scaling function map
 * \param[in]  f Signal on the sphere
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        s2let_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link s2let_analysis_px2wav
 *                        \endlink instead for complex signals.
 * \retval none
 */
void s2let_analysis_px2wav_real(
    double *f_wav,
    double *f_scal,
    const double *f,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    ssht_dl_method_t dl_method = parameters->dl_method;
    int verbosity = 0;

    complex double *flm;
    s2let_allocate_lm(&flm, L);

    switch (parameters->sampling_scheme)
    {
    case S2LET_SAMPLING_MW:
        ssht_core_mw_forward_sov_conv_sym_real(flm, f, L, dl_method, verbosity);
        break;
    case S2LET_SAMPLING_MW_SS:
        ssht_core_mw_forward_sov_conv_sym_ss_real(flm, f, L, dl_method, verbosity);
        break;
    default:
        S2LET_ERROR_GENERIC("Sampling scheme not supported.");
    }

    s2let_analysis_lm2wav_real(f_wav, f_scal, flm, parameters);

    free(flm);
}
