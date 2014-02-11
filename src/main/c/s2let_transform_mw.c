// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <complex.h>
#include <ssht.h>
#include <so3.h>
#include <stdlib.h>

/*!
 * Allocates arrays for wavelet transform in pixel space (MW sampling).
 *
 * \param[out]  f_wav Pointer to allocated space for array of wavelet
 *                    maps, using MW sampling.
 * \param[out]  f_scal Pointer to allocated space for scaling function
 *                     map, using MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \retval none
 */
void s2let_allocate_mw_f_wav(
    complex double **f_wav,
    complex double **f_scal,
    int B,
    int L,
    int J_min,
    int N
) {
    int J = s2let_j_max(L, B);
    // We actually only need N samples of the orientational angle.
    *f_wav = calloc((J-J_min+1) * (2*N-1) * L * (2*L-1), sizeof **f_wav);
    *f_scal = calloc(L * (2*L-1), sizeof **f_scal);
}

/*!
 * Spherical wavelets : full resolution analysis in real space, MW sampling.
 * Perform directional wavelet transform in real space (from scratch,
 * gives wavelet maps).
 * Sampling scheme : MW sampling.
 *
 * \param[out]  f_wav Array of wavelet maps, MW sampling.
 * \param[out]  f_scal Scaling function map, MW sampling.
 * \param[in]  f Input function (MW sampling)
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \retval none
 */
void s2let_wav_analysis_mw(
    complex double *f_wav,
    complex double *f_scal,
    const complex double *f,
    int B,
    int L,
    int J_min,
    int N
) {
    int spin = 0;
    int verbosity = 0;
    ssht_dl_method_t dl_method = SSHT_DL_TRAPANI;
    so3_parameters_t so3_parameters = {};
    so3_parameters.verbosity = verbosity;
    so3_parameters.L = L;
    so3_parameters.N = N;
    so3_parameters.n_order = SO3_N_ORDER_NEGATIVE_FIRST;
    so3_parameters.storage = SO3_STORAGE_PADDED;
    // TODO: This could probably be N_MODE_EVEN or N_MODE_ODD
    so3_parameters.n_mode = SO3_N_MODE_ALL;
    so3_parameters.dl_method = dl_method;

    int j, offset, offset_lmn;
    int J = s2let_j_max(L, B);
    //int l_min = s2let_axisym_el_min(B, J_min);

    complex double *wav_lm;
    double *scal_l;
    s2let_tiling_wavelet_allocate(&wav_lm, &scal_l, B, L, N);
    s2let_tiling_wavelet(wav_lm, scal_l, B, L, J_min, N);

    complex double *flm, *f_wav_lmn, *f_scal_lm;

    s2let_lm_allocate(&flm, L);
    ssht_core_mw_forward_sov_conv_sym(flm, f, L, spin, dl_method, verbosity);

    s2let_allocate_f_wav_lmn(&f_wav_lmn, &f_scal_lm, B, L, J_min, N);
    s2let_wav_analysis_harmonic(f_wav_lmn, f_scal_lm, flm, wav_lm, scal_l, B, L, J_min, N);

    ssht_core_mw_inverse_sov_sym(f_scal, f_scal_lm, L, spin, dl_method, verbosity);
    offset = 0;
    offset_lmn = 0;
    for (j = J_min; j <= J; ++j)
    {
        so3_core_inverse_via_ssht(
            f_wav + offset,
            f_wav_lmn + offset_lmn,
            &so3_parameters
        );
        offset_lmn += (2*N-1) * L*L;
        offset += (2*N-1) * L * (2*L-1);
    }

    free(flm);
    free(wav_lm);
    free(scal_l);
    free(f_scal_lm);
    free(f_wav_lmn);
}

/*!
 * Spherical wavelets : synthesis in real space, MW sampling.
 * Perform directional wavelet transform in real space (from scratch,
 * gives wavelet maps).
 * Sampling scheme : MW sampling.
 *
 * \param[out]  f Input function (MW sampling)
 * \param[in]  f_wav Array of wavelets maps, MW sampling.
 * \param[in]  f_scal Scaling function map, MW sampling.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \param[in]  N Azimuthal band-limit.
 * \retval none
 */
void s2let_wav_synthesis_mw(
    complex double *f,
    const complex double *f_wav,
    const complex double *f_scal,
    int B,
    int L,
    int J_min,
    int N
) {
    int spin = 0;
    int verbosity = 0;
    ssht_dl_method_t dl_method = SSHT_DL_TRAPANI;
    so3_parameters_t so3_parameters = {};
    so3_parameters.verbosity = verbosity;
    so3_parameters.L = L;
    so3_parameters.N = N;
    so3_parameters.n_order = SO3_N_ORDER_NEGATIVE_FIRST;
    so3_parameters.storage = SO3_STORAGE_PADDED;
    // TODO: This could probably be N_MODE_EVEN or N_MODE_ODD
    so3_parameters.n_mode = SO3_N_MODE_ALL;
    so3_parameters.dl_method = dl_method;

    int j, offset, offset_lmn;
    int J = s2let_j_max(L, B);
    //int l_min = s2let_axisym_el_min(B, J_min);

    complex double *wav_lm;
    double *scal_l;
    s2let_tiling_wavelet_allocate(&wav_lm, &scal_l, B, L, N);
    s2let_tiling_wavelet(wav_lm, scal_l, B, L, J_min, N);

    complex double *flm, *f_wav_lmn, *f_scal_lm;
    s2let_allocate_f_wav_lmn(&f_wav_lmn, &f_scal_lm, B, L, J_min, N);

    ssht_core_mw_forward_sov_conv_sym(f_scal_lm, f_scal, L, spin, dl_method, verbosity);
    offset = 0;
    offset_lmn = 0;
    for (j = J_min; j <= J; ++j)
    {
        so3_core_forward_via_ssht(
            f_wav_lmn + offset_lmn,
            f_wav + offset,
            &so3_parameters
        );
        offset_lmn += (2*N-1) * L*L;
        offset += (2*N-1) * L * (2*L-1);
    }

    s2let_lm_allocate(&flm, L);
    s2let_wav_synthesis_harmonic(flm, f_wav_lmn, f_scal_lm, wav_lm, scal_l, B, L, J_min, N);

    ssht_core_mw_inverse_sov_sym(f, flm, L, spin, dl_method, verbosity);

    free(flm);
    free(wav_lm);
    free(scal_l);
    free(f_scal_lm);
    free(f_wav_lmn);
}