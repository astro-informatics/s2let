// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_HPX
#define S2LET_HPX

#include <complex.h>
#include <ssht/ssht.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * Restore real healpix map from spherical harmonic coefficients.
 * Interface for HEALPIX Fortran code alm2map.
 *
 * \param[out]  f Output healpix map.
 * \param[in]  flm Spherical harmonic coefficients.
 * \param[in]  nside Healpix resolution of the output map.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_hpx_alm2map_spin_real(double* fQ, double* fU, const S2LET_COMPLEX(double)* flmE, const S2LET_COMPLEX(double)* flmB, int nside, int L, int spin);

/*!
 * Compute spherical harmonic coefficients of a real healpix map.
 * Interface for HEALPIX Fortran code map2alm.
 *
 * \param[out]  flm Spherical harmonic coefficients.
 * \param[in]  f Input healpix map.
 * \param[in]  nside Healpix resolution of the output map.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_hpx_map2alm_spin_real(S2LET_COMPLEX(double)* flmE, S2LET_COMPLEX(double)* flmB, const double* fQ, const double* fU, int nside, int L, int spin);

/*!
 * Restore real healpix map from spherical harmonic coefficients.
 * Interface for HEALPIX Fortran code alm2map.
 *
 * \param[out]  f Output healpix map.
 * \param[in]  flm Spherical harmonic coefficients.
 * \param[in]  nside Healpix resolution of the output map.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_hpx_alm2map_real(double* f, const S2LET_COMPLEX(double)* flm, int nside, int L);


/*!
 * Compute spherical harmonic coefficients of a real healpix map.
 * Interface for HEALPIX Fortran code map2alm.
 *
 * \param[out]  flm Spherical harmonic coefficients.
 * \param[in]  f Input healpix map.
 * \param[in]  nside Healpix resolution of the output map.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
void s2let_hpx_map2alm_real(S2LET_COMPLEX(double)* flm, const double* f, int nside, int L);


/*!
 * Read Healpix map from a FITS file.
 * Interface for HEALPIX Fortran code read_bintab.
 *
 * \param[out]  f Input healpix map.
 * \param[in]  file Filename.
 * \param[in]  nside Healpix resolution of the output map.
 * \retval none
 */
void s2let_hpx_read_map(double* f, char* file, int nside);

void s2let_hpx_read_maps(double* f, char* file, int nside, int nmaps);


/*!
 * Write Healpix map to a FITS file.
 * Interface for HEALPIX Fortran code write_bintab.
 *
 * \param[in]  file Filename.
 * \param[in]  f Input healpix map.
 * \param[in]  nside Healpix resolution of the output map.
 * \retval none
 */
void s2let_hpx_write_map(char* file, const double* f, int nside);


/*!
 * Allocate Healpix map.
 *
 * \param[inout]  f healpix map.
 * \param[in]  nside Healpix resolution.
 * \retval none
 */
void s2let_hpx_allocate_real(double **f, int nside);

#ifdef __cplusplus
}
#endif
#endif
