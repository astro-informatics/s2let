// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_AXISYM_HPX
#define S2LET_AXISYM_HPX


#ifdef __cplusplus
extern "C" {
#endif


/*!
 * Allocates arrays for wavelets and scaling functions in pixel space (Healpix sampling).
 *
 * \param[out]  f_wav Array of wavelets HEALPIX maps.
 * \param[out]  f_scal Scaling function HEALPIX map.
 * \param[in]  nside HEALPIX resolution parameter.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_allocate_hpx_f_wav_real(double **f_wav, double **f_scal, int nside, const s2let_parameters_t *parameters);

/*!
 * Perform wavelet transform in real space (from scratch, gives wavelet maps).
 * Input function is real.
 * Sampling scheme : HEALPIX sampling.
 * Spherical wavelets : analysis in real space, HEALPIX sampling.
 * Note : multiresolution in used in harmonic space but all maps are at resolution nside.
 *
 * \param[out]  f_wav Array of wavelets HEALPIX maps.
 * \param[out]  f_scal Scaling function HEALPIX map.
 * \param[in]  f Input function (HEALPIX map)
 * \param[in]  nside HEALPIX resolution parameter.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_wav_analysis_hpx_real(double *f_wav, double *f_scal, const double *f, int nside, const s2let_parameters_t *parameters);

/*!
 * Perform wavelet transform in real space (from scratch, gives wavelet maps).
 * Input function is real.
 * Sampling scheme : HEALPIX sampling.
 * Spherical wavelets : synthesis in real space, HEALPIX sampling.
 * Note : multiresolution in used in harmonic space but all maps are at resolution nside.
 *
 * \param[out]  f Input function (HEALPIX map)
 * \param[in]  f_wav Array of wavelets HEALPIX maps.
 * \param[in]  f_scal Scaling function HEALPIX map.
 * \param[in]  nside HEALPIX resolution parameter.
 * \param[in]  B Wavelet parameter.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  J_min First wavelet scale to be used.
 * \retval none
 */
void s2let_transform_axisym_wav_synthesis_hpx_real(double *f, const double *f_wav, const double *f_scal, int nside, const s2let_parameters_t *parameters);

#ifdef __cplusplus
}
#endif
#endif