// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_ANALYSIS_ADJOINT
#define S2LET_ANALYSIS_ADJOINT

#include <ssht/ssht.h>

#ifdef __cplusplus
extern "C" {
#endif


//! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*!
 * Wavelet analysis adjoint from Wigner space to harmonic space for complex signals.
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
void s2let_analysis_adjoint_lmn2lm(
    S2LET_COMPLEX(double) *flm,
    const S2LET_COMPLEX(double) *f_wav_lmn,
    const S2LET_COMPLEX(double) *f_scal_lm,
    const S2LET_COMPLEX(double) *wav_lm,
    const double *scal_l,
    const s2let_parameters_t *parameters);



//! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


/*!
 * Wavelet analysis adjoint from Wigner space to harmonic space for real signals.
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
void s2let_analysis_adjoint_lmn2lm_real(
    S2LET_COMPLEX(double) *flm,
    const S2LET_COMPLEX(double) *f_wav_lmn,
    const S2LET_COMPLEX(double) *f_scal_lm,
    const S2LET_COMPLEX(double) *wav_lm,
    const double *scal_l,
    const s2let_parameters_t *parameters);




//! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*!
 * Wavelet analysis adjoint from wavelet space to harmonic space for complex signals.
 * with fully manual wavelet tiling, using multiresolution as default with 
 * the band-limits provided in input.
 *
 * \param[out]  flm Spherical harmonic coefficients of the signal
 * \param[in]  f_wav Array of wavelets maps
 * \param[in]  f_scal Scaling function map
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
void s2let_analysis_adjoint_wav2lm_manual(
    S2LET_COMPLEX(double) *flm,
    const S2LET_COMPLEX(double) *f_wav,
    const S2LET_COMPLEX(double) *f_scal,
    const double *scal_l,
    const S2LET_COMPLEX(double) *wav_lm,
    const int scal_bandlimit,
    const int *wav_bandlimits,
    int J,
    int L,
    int spin,
    int N );



//! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*!
 * Wavelet analysis adjoint from wavelet space to harmonic space for complex signals.
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
void s2let_analysis_adjoint_wav2lm(
    S2LET_COMPLEX(double) *flm,
    const S2LET_COMPLEX(double) *f_wav,
    const S2LET_COMPLEX(double) *f_scal,
    const s2let_parameters_t *parameters );



//! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*!
 * Wavelet analysis adjoint from wavelet space to harmonic space for real signals.
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
void s2let_analysis_adjoint_wav2lm_real(
    S2LET_COMPLEX(double) *flm,
    const double *f_wav,
    const double *f_scal,
    const s2let_parameters_t *parameters);


//! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*!
 * Wavelet analysis adjoint from wavelet space to pixel space for complex signals.
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
void s2let_analysis_adjoint_wav2px(
    S2LET_COMPLEX(double) *f,
    const S2LET_COMPLEX(double) *f_wav,
    const S2LET_COMPLEX(double) *f_scal,
    const s2let_parameters_t *parameters );

//! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
void s2let_analysis_adjoint_wav2px_real(
    double *f,
    const double *f_wav,
    const double *f_scal,
    const s2let_parameters_t *parameters );





#ifdef __cplusplus
}
#endif
#endif
