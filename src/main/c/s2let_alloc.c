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
 * Allocate map for a complex signal in pixel space using MW sampling.
 * \param[out]  f Allocated map
 * \param[in]   L Harmonic band-limit
 */
void s2let_allocate_mw(complex double **f, int L)
{
  *f = calloc(L * (2*L-1), sizeof **f);
}

/*!
 * Allocate map for a real signal in pixel space using MW sampling.
 * \param[out]  f Allocated map
 * \param[in]   L Harmonic band-limit
 */
void s2let_allocate_mw_real(double **f, int L)
{
  *f = calloc(L * (2*L-1), sizeof **f);
}

/*!
 * Allocate map for a complex signal in pixel space using MWSS sampling.
 * \param[out]  f Allocated map
 * \param[in]   L Harmonic band-limit
 */
void s2let_allocate_mwss(complex double **f, int L)
{
  *f = calloc((2*L)*(L+1), sizeof **f);
}

/*!
 * Allocate map for a real signal in pixel space using MWSS sampling.
 * \param[out]  f Allocated map
 * \param[in]   L Harmonic band-limit
 */
void s2let_allocate_mwss_real(double **f, int L)
{
  *f = calloc((2*L)*(L+1), sizeof **f);
}

/*!
 * Allocate spherical harmonic coefficients for a given
 * bandlimit L.
 *
 * \param[out]  flm Pointer to allocated space for spherical
 *                  harmonic coefficients.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_allocate_lm(complex double **flm, int L)
{
    *flm = calloc(L * L, sizeof **flm);
}



/*!
 * Allocates arrays for directional wavelet transform in Wigner space.
 *
 * \param[out]  f_wav_lmn Wigner coefficients of the wavelet contributions.
 *                        Each wavelet has size (2*N-1)*L*L and there are
 *                        (J-J_min+1) scales.
 * \param[out]  f_scal_lm Spherical harmonic coefficients of the scaling
 *                        contribution (L*L).
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link s2let_parameters_t::B B\endlink,
 *                        \link s2let_parameters_t::L L\endlink,
 *                        \link s2let_parameters_t::J_min J_min\endlink
 *                        \link s2let_parameters_t::N N\endlink
 * \retval none
 */
void s2let_allocate_lmn_f_wav(
    complex double **f_wav_lmn,
    complex double **f_scal_lm,
    const s2let_parameters_t *parameters
) {
    *f_wav_lmn = calloc(s2let_n_lmn_wav(parameters), sizeof **f_wav_lmn);
    *f_scal_lm = calloc(s2let_n_lm_scal(parameters), sizeof **f_scal_lm);
}

/*!
 * Allocates arrays for wavelet transform in wavelet space.
 *
 * \param[out]  f_wav Pointer to allocated space for array of wavelet
 *                    maps, using MW sampling.
 * \param[out]  f_scal Pointer to allocated space for scaling function
 *                     map, using MW sampling.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link s2let_parameters_t::B B\endlink,
 *                        \link s2let_parameters_t::L L\endlink,
 *                        \link s2let_parameters_t::J_min J_min\endlink
 *                        \link s2let_parameters_t::N N\endlink
 * \retval none
 */
void s2let_allocate_f_wav(
    complex double **f_wav,
    complex double **f_scal,
    const s2let_parameters_t *parameters
) {
    *f_wav = calloc(s2let_n_wav(parameters), sizeof **f_wav);
    *f_scal = calloc(s2let_n_scal(parameters), sizeof **f_scal);
}

/*!
 * Allocates arrays for wavelet transform of real signal in wavelet space.
 *
 * \param[out]  f_wav Pointer to allocated space for array of wavelet
 *                    maps, using MW sampling.
 * \param[out]  f_scal Pointer to allocated space for scaling function
 *                     map, using MW sampling.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link s2let_parameters_t::B B\endlink,
 *                        \link s2let_parameters_t::L L\endlink,
 *                        \link s2let_parameters_t::J_min J_min\endlink
 *                        \link s2let_parameters_t::N N\endlink
 * \retval none
 */
void s2let_allocate_f_wav_real(
    double **f_wav,
    double **f_scal,
    const s2let_parameters_t *parameters
) {
    *f_wav = calloc(s2let_n_wav(parameters), sizeof **f_wav);
    *f_scal = calloc(s2let_n_scal(parameters), sizeof **f_scal);
}




// CURVELETS:
/*!
 * Allocates arrays for curvelet transform in Wigner space.
 *
 * \param[out]  f_cur_lmn Wigner coefficients of the curvelet contributions.
 *                        Each wavelet has size (2*N-1)*L*L and there are
 *                        (J-J_min+1) scales.
 * \param[out]  f_scal_lm Spherical harmonic coefficients of the scaling
 *                        contribution (L*L).
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link s2let_parameters_t::B B\endlink,
 *                        \link s2let_parameters_t::L L\endlink,
 *                        \link s2let_parameters_t::J_min J_min\endlink
 *                        \link s2let_parameters_t::N N\endlink
 * \retval none
 */
void s2let_allocate_lmn_f_cur(
    complex double **f_cur_lmn,
    complex double **f_scal_lm,
    const s2let_parameters_t *parameters
) {
    *f_cur_lmn = calloc(s2let_n_lmn_cur(parameters), sizeof **f_cur_lmn);
    *f_scal_lm = calloc(s2let_n_lm_scal(parameters), sizeof **f_scal_lm);
    
}
/*!
 * Allocates arrays for curvelet transform in curvelet space.
 *
 * \param[out]  f_cur Pointer to allocated space for array of curvelet
 *                    maps, using MW sampling.
 * \param[out]  f_scal Pointer to allocated space for scaling function
 *                     map, using MW sampling.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link s2let_parameters_t::B B\endlink,
 *                        \link s2let_parameters_t::L L\endlink,
 *                        \link s2let_parameters_t::J_min J_min\endlink
 *                        \link s2let_parameters_t::N N\endlink
 * \retval none
 */
void s2let_allocate_f_cur(
     complex double **f_cur,
     complex double **f_scal,
     const s2let_parameters_t *parameters
) {
    *f_cur = calloc(s2let_n_cur(parameters), sizeof **f_cur);
    *f_scal = calloc(s2let_n_scal(parameters), sizeof **f_scal);
}

/*!
 * Allocates arrays for wavelet transform of real signal in curvelet space.
 *
 * \param[out]  f_cur Pointer to allocated space for array of curvelet
 *                    maps, using MW sampling.
 * \param[out]  f_scal Pointer to allocated space for scaling function
 *                     map, using MW sampling.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link s2let_parameters_t::B B\endlink,
 *                        \link s2let_parameters_t::L L\endlink,
 *                        \link s2let_parameters_t::J_min J_min\endlink
 *                        \link s2let_parameters_t::N N\endlink
 * \retval none
 */
void s2let_allocate_f_cur_real(
    double **f_cur,
    double **f_scal,
    const s2let_parameters_t *parameters
) {
    *f_cur = calloc(s2let_n_cur(parameters), sizeof **f_cur);
    *f_scal = calloc(s2let_n_scal(parameters), sizeof **f_scal);
}
