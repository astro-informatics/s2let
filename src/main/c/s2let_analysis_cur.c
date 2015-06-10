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
 * Curvelet analysis from harmonic space to Wigner space for complex signals.
 *
 * \param[out]  f_cur_lmn Curvelet transform (Wigner coefficients of curvelet contribution).
 * \param[out]  f_scal_lm Curvelet transform (Spherical harmonic coefficients of scaling contribution).
 * \param[in]  flm Spherical harmonic coefficients of input function.
 * \param[in]  cur_lm Curvelet kernels in harmonic space.
 * \param[in]  scal_l Scaling function kernels in harmonic space.
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        s2let_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link s2let_analysis_cur_lm2lmn_real
 *                        \endlink instead for real signals.
 * \retval none
 */
void s2let_analysis_cur_lm2lmn(
    complex double *f_cur_lmn,
    complex double *f_scal_lm,
    const complex double *flm,
    const complex double *cur_lm,
    const double *scal_l,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = L ;
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
    lmn_ind = 0;
    lm_ind = 0;
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

        for (n = -Nj+1; n < Nj; n+=1)
        {
            for (el = MAX(ABS(spin), ABS(n)); el < bandlimit; ++el)
            {
                ssht_sampling_elm2ind(&lm_ind, el, n);
                psi = 8.*PI*PI/(2.*el+1.) * conj(cur_lm[j*L*L + lm_ind]);
                for (m = -el; m <= el; ++m)
                {
                    ssht_sampling_elm2ind(&lm_ind, el, m);
                    so3_sampling_elmn2ind(&lmn_ind, el, m, n, &so3_parameters);
                    f_cur_lmn[offset + lmn_ind] = flm[lm_ind] * psi;
                    
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
 * Curvelet analysis from harmonic space to Wigner space for real signals.
 *
 * \param[out]  f_cur_lmn Curvelet transform (Wigner coefficients of curvelet contribution).
 * \param[out]  f_scal_lm Curvelet transform (spherical harmonic coefficients of scaling contribution).
 * \param[in]  flm Spherical harmonic coefficients of input function.
 * \param[in]  cur_lm Curvelet kernels.
 * \param[in]  scal_l Scaling function kernels.
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        s2let_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link s2let_analysis_cur_lm2lmn
 *                        \endlink instead for complex signals.
 * \retval none
 */
void s2let_analysis_cur_lm2lmn_real(
    complex double *f_cur_lmn,
    complex double *f_scal_lm,
    const complex double *flm,
    const complex double *cur_lm,
    const double *scal_l,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = L;  //parameters->N;

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

        for (n = 1-Nj%2; n < Nj; n+=1)
        {
            for (el = n; el < bandlimit; ++el)
            {
                ssht_sampling_elm2ind(&lm_ind, el, n);
                psi = 8*PI*PI/(2*el+1) * conj(cur_lm[j*L*L + lm_ind]);
                for (m = -el; m <= el; ++m)
                {
                    ssht_sampling_elm2ind(&lm_ind, el, m);
                    so3_sampling_elmn2ind_real(&lmn_ind, el, m, n, &so3_parameters);
                    f_cur_lmn[offset + lmn_ind] = flm[lm_ind] * psi;
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
 * curvelet analysis from harmonic space to curvelet space for complex signals.
 *
 * \param[out]  f_cur Array of curvelet maps
 * \param[out]  f_scal Scaling function map
 * \param[in]  flm Spherical harmonic coefficients of the signal
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        s2let_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link s2let_analysis_lm2cur_real
 *                        \endlink instead for real signals.
 * \retval none
 */
void s2let_analysis_lm2cur(
                           complex double *f_cur,
                           complex double *f_scal,
                           const complex double *flm,
                           const s2let_parameters_t *parameters
                           ) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = parameters->N;  //L;
    ssht_dl_method_t dl_method = parameters->dl_method;
    
    int bandlimit = L;
    int verbosity = 0;
    so3_parameters_t so3_parameters = {};
    fill_so3_parameters(&so3_parameters, parameters);
    
    int j, offset_cur, offset_cur_lmn;
    int J = s2let_j_max(parameters);
    
    complex double *cur_lm;
    double *scal_l;
    s2let_tiling_curvelet_allocate(&cur_lm, &scal_l, parameters);
    s2let_tiling_curvelet(cur_lm, scal_l, parameters);
    
    complex double *f_cur_lmn, *f_scal_lm;
    
    s2let_allocate_lmn_f_cur(&f_cur_lmn, &f_scal_lm, parameters);
    s2let_analysis_cur_lm2lmn(f_cur_lmn, f_scal_lm, flm, cur_lm, scal_l, parameters);
    
    
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
    
    offset_cur = 0;
    offset_cur_lmn = 0;
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
        
        so3_core_inverse_via_ssht(f_cur + offset_cur,
                                  f_cur_lmn + offset_cur_lmn,
                                  &so3_parameters
                                  );

        
        offset_cur_lmn += so3_sampling_flmn_size(&so3_parameters);
        offset_cur += so3_sampling_f_size(&so3_parameters);
        
    }
    
    
    free(cur_lm);
    free(scal_l);
    free(f_scal_lm);
    free(f_cur_lmn);
}



/*!
 * Curvelet analysis from harmonic space to curvelet space for real signals.
 *
 * \param[out]  f_cur Array of curvelet maps
 * \param[out]  f_scal Scaling function map
 * \param[in]  flm Spherical harmonic coefficients of the signal
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        s2let_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link s2let_analysis_lm2cur
 *                        \endlink instead for complex signals.
 * \retval none
 */
void s2let_analysis_lm2cur_real(
    double *f_cur,
    double *f_scal,
    const complex double *flm,
    const s2let_parameters_t *parameters
) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int N = L ; //parameters->N;
    ssht_dl_method_t dl_method = parameters->dl_method;

    s2let_parameters_t real_parameters = *parameters;
    real_parameters.reality = 1;

    int bandlimit = L;
    int verbosity = 0;
    so3_parameters_t so3_parameters = {};
    fill_so3_parameters(&so3_parameters, &real_parameters);

    int j, offset, offset_lmn;
    int J = s2let_j_max(&real_parameters);

    complex double *cur_lm;
    double *scal_l;
    s2let_tiling_curvelet_allocate(&cur_lm, &scal_l, &real_parameters);
    s2let_tiling_curvelet(cur_lm, scal_l, &real_parameters);

    complex double *f_cur_lmn, *f_scal_lm;

    s2let_allocate_lmn_f_cur(&f_cur_lmn, &f_scal_lm, &real_parameters);
    s2let_analysis_cur_lm2lmn_real(f_cur_lmn, f_scal_lm, flm, cur_lm, scal_l, &real_parameters);

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
            f_cur + offset,
            f_cur_lmn + offset_lmn,
            &so3_parameters
        );
        offset_lmn += so3_sampling_flmn_size(&so3_parameters);
        offset += so3_sampling_f_size(&so3_parameters);
    }

    free(cur_lm);
    free(scal_l);
    free(f_scal_lm);
    free(f_cur_lmn);
}

/*!
 * Curvelet analysis from pixel space to curvelet space for complex signals.
 *
 * \param[out]  f_cur Array of curvelet maps
 * \param[out]  f_scal Scaling function map
 * \param[in]  f Signal on the sphere
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        s2let_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link s2let_analysis_px2cur_real
 *                        \endlink instead for real signals.
 * \retval none
 */
void s2let_analysis_px2cur(
    complex double *f_cur,
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

    s2let_analysis_lm2cur(f_cur, f_scal, flm, parameters);

    free(flm);
}

/*!
 * Curvelet analysis from pixel space to curvelet space for real signals.
 *
 * \param[out]  f_cur Array of curvelet maps
 * \param[out]  f_scal Scaling function map
 * \param[in]  f Signal on the sphere
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        s2let_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link s2let_analysis_px2cur
 *                        \endlink instead for complex signals.
 * \retval none
 */
void s2let_analysis_px2cur_real(
    double *f_cur,
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

    s2let_analysis_lm2cur_real(f_cur, f_scal, flm, parameters);

    free(flm);
}