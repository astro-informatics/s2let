// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_FITS
#define S2LET_FITS

/*!
 * Read Healpix resolution from a FITS file.
 *
 * \param[in]  file Filename.
 * \retval int resolution parameter
 */
int s2let_fits_hpx_read_nside(char* filename);

/*!
 * Read MW resolution / band-limit parameter from a FITS file.
 *
 * \param[in]  file Filename.
 * \retval int resolution parameter
 */
int s2let_fits_mw_read_bandlimit(char* filename);

/*!
 * Read MW map from a FITS file.
 *
 * \param[out]  f Input map (MW sampling).
 * \param[in]  file Filename.
 * \param[in]  L Band-limit / resolution parameter.
 * \retval none
 */
void s2let_fits_mw_read_map(double* f, char* file, int L);

/*!
 * Read MWSS map from a FITS file.
 *
 * \param[out]  f Input map (MW sampling).
 * \param[in]  file Filename.
 * \param[in]  L Band-limit / resolution parameter.
 * \retval none
 */
void s2let_fits_mwss_read_map(double* f, char* file, int L);

/*!
 * Write MW map from a FITS file.
 *
 * \param[in]  f Input map (MW sampling).
 * \param[in]  file Filename.
 * \param[in]  L Band-limit / resolution parameter.
 * \retval none
 */
void s2let_fits_mw_write_map(char* file, double* f, int L);

/*!
 * Write MWSS map from a FITS file.
 *
 * \param[in]  f Input map (MW sampling).
 * \param[in]  file Filename.
 * \param[in]  L Band-limit / resolution parameter.
 * \retval none
 */
void s2let_fits_mwss_write_map(char* file, double* f, int L);

/*!
 * Read MW map from a FITS file.
 *
 * \param[out]  f Input map (MW sampling).
 * \param[in]  file Filename.
 * \param[in]  L Band-limit / resolution parameter.
 * \retval none
 */
void s2let_fits_mw_read_spin_maps(double* fQ, double* fU, char* file, int L);

/*!
 * Read MWSS map from a FITS file.
 *
 * \param[out]  f Input map (MW sampling).
 * \param[in]  file Filename.
 * \param[in]  L Band-limit / resolution parameter.
 * \retval none
 */
void s2let_fits_mwss_read_spin_maps(double* fQ, double* fU, char* file, int L);

/*!
 * Write MW map from a FITS file.
 *
 * \param[in]  f Input map (MW sampling).
 * \param[in]  file Filename.
 * \param[in]  L Band-limit / resolution parameter.
 * \retval none
 */
void s2let_fits_mw_write_spin_maps(char* file, double* fQ, double* fU, int L);

/*!
 * Write MWSS map from a FITS file.
 *
 * \param[in]  f Input map (MW sampling).
 * \param[in]  file Filename.
 * \param[in]  L Band-limit / resolution parameter.
 * \retval none
 */
void s2let_fits_mwss_write_spin_maps(char* file, double* fQ, double* fU, int L);

#endif
