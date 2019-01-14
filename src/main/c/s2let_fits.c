// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include "fitsio.h"
#include "s2let.h"


void printerror(int status)
{
  if (status){
    fits_report_error(stderr, status); /* print error report */
    exit( status );    /* terminate the program, returning error status */
  }
  return;
}

int s2let_fits_mw_read_bandlimit(char* filename)
{
  long     naxes, *naxis, Lread;
  int      status, hdutype, nfound;
  char     comment[FLEN_COMMENT];
  fitsfile *fptr;
  status = 0;


  if ( fits_open_file(&fptr, filename, READONLY, &status) )
    printerror( status );

  if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) )
    printerror( status );

  if (hdutype != BINARY_TBL)
    fprintf(stderr, "%s (%d): Extension is not binary!\n", __FILE__, __LINE__);

  if ( fits_read_key_lng(fptr, "NAXIS", &naxes, comment, &status) )
    printerror( status );

  naxis = (long *)malloc(((size_t)naxes)*sizeof(long));
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, naxes, naxis, &nfound, &status)
       || nfound != naxes )
    printerror( status );

  if ( fits_read_key_lng(fptr, "L", &Lread, comment, &status) )
  {
    Lread = -1;
    status = 0;
  }


  if ( fits_close_file(fptr, &status) )
    printerror( status );

  return Lread;

}


int s2let_fits_hpx_read_nside(char* filename)
{
  long     naxes, *naxis, nsideread;
  int      status, hdutype, nfound;
  char     comment[FLEN_COMMENT];
  fitsfile *fptr;
  status = 0;


  if ( fits_open_file(&fptr, filename, READONLY, &status) )
    printerror( status );

  if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) )
    printerror( status );

  if (hdutype != BINARY_TBL)
    fprintf(stderr, "%s (%d): Extension is not binary!\n", __FILE__, __LINE__);

  if ( fits_read_key_lng(fptr, "NAXIS", &naxes, comment, &status) )
    printerror( status );

  naxis = (long *)malloc(((size_t)naxes)*sizeof(long));
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, naxes, naxis, &nfound, &status)
       || nfound != naxes )
    printerror( status );

  if ( fits_read_key_lng(fptr, "NSIDE", &nsideread, comment, &status) )
  {
    nsideread = -1;
    status = 0;
  }

  if ( fits_close_file(fptr, &status) )
    printerror( status );

  return nsideread;

}

void s2let_fits_mw_write_map(char* filename, double* f, int L)
{
  fitsfile *fptr;
  int status, hdutype;
  long firstrow, firstelem;
  int bitpix   =  SHORT_IMG;
  long naxis   =   0;
  long naxes[] = {0,0};
  int tfields   = 1;
  char extname[] = "BINTABLE";
  char *ttype[] = { "SIGNAL" };
  char *tform[] = { "1D" };
  char *tunit[] = { " " };

  long npix = L * (2*L-1);
  status = 0;

  if (fits_create_file(&fptr, filename, &status))
    fprintf(stderr, "%s (%d): Could not create new fits file.\n",
      __FILE__, __LINE__);

  if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
    fprintf(stderr, "%s (%d): Could not create new image file.\n",
      __FILE__, __LINE__);

  if ( fits_movabs_hdu(fptr, 1, &hdutype, &status) )
    fprintf(stderr, "%s (%d): Could not move to first HDU.\n",
      __FILE__, __LINE__);

  if ( fits_create_tbl( fptr, BINARY_TBL, npix, tfields, ttype, tform,
      tunit, extname, &status) )
    fprintf(stderr, "%s (%d): Could not create new binary table.\n",
      __FILE__, __LINE__);

  if (fits_write_key(fptr, TINT, "L", &L,
         "Resolution parameter", &status))
    fprintf(stderr, "%s (%d): Could not write L keyword.\n",
      __FILE__, __LINE__);

  firstrow  = 1;
  firstelem = 1;

  if (fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, npix, f, &status))
    fprintf(stderr, "%s (%d): Could not write signal.\n", __FILE__, __LINE__);

  if ( fits_close_file(fptr, &status) )       /* close the FITS file */
    fprintf(stderr, "%s (%d): Could not close file.\n",
      __FILE__, __LINE__);

}

void s2let_fits_mwss_write_map(char* filename, double* f, int L)
{
  fitsfile *fptr;
  int status, hdutype;
  long firstrow, firstelem;
  int bitpix   =  SHORT_IMG;
  long naxis   =   0;
  long naxes[] = {0,0};
  int tfields   = 1;
  char extname[] = "BINTABLE";
  char *ttype[] = { "SIGNAL" };
  char *tform[] = { "1D" };
  char *tunit[] = { " " };

  long npix = (L+1) * 2*L;
  status = 0;

  if (fits_create_file(&fptr, filename, &status))
    fprintf(stderr, "%s (%d): Could not create new fits file.\n",
      __FILE__, __LINE__);

  if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
    fprintf(stderr, "%s (%d): Could not create new image file.\n",
      __FILE__, __LINE__);

  if ( fits_movabs_hdu(fptr, 1, &hdutype, &status) )
    fprintf(stderr, "%s (%d): Could not move to first HDU.\n",
      __FILE__, __LINE__);

  if ( fits_create_tbl( fptr, BINARY_TBL, npix, tfields, ttype, tform,
      tunit, extname, &status) )
    fprintf(stderr, "%s (%d): Could not create new binary table.\n",
      __FILE__, __LINE__);

  if (fits_write_key(fptr, TINT, "L", &L,
         "Resolution parameter", &status))
    fprintf(stderr, "%s (%d): Could not write L keyword.\n",
      __FILE__, __LINE__);

  int mwss = 1;
  if (fits_write_key(fptr, TINT, "MWSS", &mwss,
         "Sampling scheme parameter", &status))
    fprintf(stderr, "%s (%d): Could not write MWSS keyword.\n",
      __FILE__, __LINE__);

  firstrow  = 1;
  firstelem = 1;

  if (fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, npix, f, &status))
    fprintf(stderr, "%s (%d): Could not write signal.\n", __FILE__, __LINE__);

  if ( fits_close_file(fptr, &status) )       /* close the FITS file */
    fprintf(stderr, "%s (%d): Could not close file.\n",
      __FILE__, __LINE__);

}


void s2let_fits_mw_read_map(double* f, char* filename, int L)
{
  long     naxes, *naxis, npix, npercol, irow, Lread;
  int      status, hdutype, nfound, anynul;
  float    nulval;
  char     comment[FLEN_COMMENT];
  fitsfile *fptr;
  status = 0;

  if ( fits_open_file(&fptr, filename, READONLY, &status) )
    printerror( status );

  if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) )
    printerror( status );


  if (hdutype != BINARY_TBL)
    fprintf(stderr, "%s (%d): Extension is not binary!\n", __FILE__, __LINE__);

  if ( fits_read_key_lng(fptr, "NAXIS", &naxes, comment, &status) )
    printerror( status );

  naxis = (long *)malloc(((size_t)naxes)*sizeof(long));
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, naxes, naxis, &nfound, &status)
       || nfound != naxes )
    printerror( status );

  if ( fits_read_key_lng(fptr, "L", &Lread, comment, &status) )
    printerror(status);

  if( Lread != L )
    printf("Attention : read L = %li but you specified L = %i\n",Lread,L);

  npix = L * (2*L-1);

  if ( (npix%naxis[1]) != 0 )
    fprintf(stderr, "%s (%d): Problem with npix.\n", __FILE__, __LINE__);

  npercol = npix/naxis[1];
  nulval = -1.6375e30;
  for (irow = 0; irow < naxis[1]; irow++) {

    if ( fits_read_col(fptr, TDOUBLE, 1, irow+1, 1, npercol, &nulval,
		       &(f[irow*npercol]), &anynul, &status) ) {
      printerror(status);
    }
  }

  if ( fits_close_file(fptr, &status) )
    printerror( status );

}

void s2let_fits_mwss_read_map(double* f, char* filename, int L)
{
  long     naxes, *naxis, npix, npercol, irow, Lread;
  int      status, hdutype, nfound, anynul;
  float    nulval;
  char     comment[FLEN_COMMENT];
  fitsfile *fptr;
  status = 0;

  if ( fits_open_file(&fptr, filename, READONLY, &status) )
    printerror( status );

  if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) )
    printerror( status );


  if (hdutype != BINARY_TBL)
    fprintf(stderr, "%s (%d): Extension is not binary!\n", __FILE__, __LINE__);

  if ( fits_read_key_lng(fptr, "NAXIS", &naxes, comment, &status) )
    printerror( status );

  naxis = (long *)malloc(((size_t)naxes)*sizeof(long));
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, naxes, naxis, &nfound, &status)
       || nfound != naxes )
    printerror( status );

  if ( fits_read_key_lng(fptr, "L", &Lread, comment, &status) )
    printerror(status);

  if( Lread != L )
    printf("Attention : read L = %li but you specified L = %i\n",Lread,L);


  long mwss_read;
  if ( fits_read_key_lng(fptr, "MWSS", &mwss_read, comment, &status) )
    printerror(status);

  if( mwss_read != 1 )
    printf("Attention : the file is not stored with MWSS sampling.\n");

  npix = (L+1) * 2*L;

  if ( (npix%naxis[1]) != 0 )
    fprintf(stderr, "%s (%d): Problem with npix.\n", __FILE__, __LINE__);

  npercol = npix/naxis[1];
  nulval = -1.6375e30;
  for (irow = 0; irow < naxis[1]; irow++) {

    if ( fits_read_col(fptr, TDOUBLE, 1, irow+1, 1, npercol, &nulval,
               &(f[irow*npercol]), &anynul, &status) ) {
      printerror(status);
    }
  }

  if ( fits_close_file(fptr, &status) )
    printerror( status );

}


void s2let_fits_mw_write_spin_maps(char* filename, double* fQ, double*fU, int L)
{
  fitsfile *fptr;
  int status, hdutype;
  long firstrow, firstelem;
  int bitpix   =  SHORT_IMG;
  long naxis   =   0;
  long naxes[] = {0,0};
  int tfields   = 2;
  char extname[] = "BINTABLE";
  char *ttype[] = { "Q_POLARISATION", "U_POLARISATION" };
  char *tform[] = { "1D", "1D" };
  char *tunit[] = { " ", " " };

  long npix = L * (2*L-1);
  status = 0;

  if (fits_create_file(&fptr, filename, &status))
    fprintf(stderr, "%s (%d): Could not create new fits file.\n",
	    __FILE__, __LINE__);

  if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
    fprintf(stderr, "%s (%d): Could not create new image file.\n",
	    __FILE__, __LINE__);

  if ( fits_movabs_hdu(fptr, 1, &hdutype, &status) )
    fprintf(stderr, "%s (%d): Could not move to first HDU.\n",
	    __FILE__, __LINE__);

  if ( fits_create_tbl( fptr, BINARY_TBL, npix, tfields, ttype, tform,
			tunit, extname, &status) )
    fprintf(stderr, "%s (%d): Could not create new binary table.\n",
	    __FILE__, __LINE__);

  if (fits_write_key(fptr, TINT, "L", &L,
		     "Resolution parameter", &status))
    fprintf(stderr, "%s (%d): Could not write L keyword.\n",
	    __FILE__, __LINE__);

  firstrow  = 1;
  firstelem = 1;

  if (fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, npix, fQ, &status))
    fprintf(stderr, "%s (%d): Could not write signal.\n", __FILE__, __LINE__);

  if (fits_write_col(fptr, TDOUBLE, 2, firstrow, firstelem, npix, fU, &status))
    fprintf(stderr, "%s (%d): Could not write signal.\n", __FILE__, __LINE__);

  if ( fits_close_file(fptr, &status) )       /* close the FITS file */
    fprintf(stderr, "%s (%d): Could not close file.\n",
	    __FILE__, __LINE__);

}

void s2let_fits_mwss_write_spin_maps(char* filename, double* fQ, double*fU, int L)
{
  fitsfile *fptr;
  int status, hdutype;
  long firstrow, firstelem;
  int bitpix   =  SHORT_IMG;
  long naxis   =   0;
  long naxes[] = {0,0};
  int tfields   = 2;
  char extname[] = "BINTABLE";
  char *ttype[] = { "Q_POLARISATION", "U_POLARISATION" };
  char *tform[] = { "1D", "1D" };
  char *tunit[] = { " ", " " };

  long npix = (L+1) * 2*L;
  status = 0;

  if (fits_create_file(&fptr, filename, &status))
    fprintf(stderr, "%s (%d): Could not create new fits file.\n",
        __FILE__, __LINE__);

  if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
    fprintf(stderr, "%s (%d): Could not create new image file.\n",
        __FILE__, __LINE__);

  if ( fits_movabs_hdu(fptr, 1, &hdutype, &status) )
    fprintf(stderr, "%s (%d): Could not move to first HDU.\n",
        __FILE__, __LINE__);

  if ( fits_create_tbl( fptr, BINARY_TBL, npix, tfields, ttype, tform,
            tunit, extname, &status) )
    fprintf(stderr, "%s (%d): Could not create new binary table.\n",
        __FILE__, __LINE__);

  if (fits_write_key(fptr, TINT, "L", &L,
             "Resolution parameter", &status))
    fprintf(stderr, "%s (%d): Could not write L keyword.\n",
        __FILE__, __LINE__);

  int mwss = 1;
  if (fits_write_key(fptr, TINT, "MWSS", &mwss,
         "Sampling scheme parameter", &status))
    fprintf(stderr, "%s (%d): Could not write MWSS keyword.\n",
      __FILE__, __LINE__);

  firstrow  = 1;
  firstelem = 1;

  if (fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, npix, fQ, &status))
    fprintf(stderr, "%s (%d): Could not write signal.\n", __FILE__, __LINE__);

  if (fits_write_col(fptr, TDOUBLE, 2, firstrow, firstelem, npix, fU, &status))
    fprintf(stderr, "%s (%d): Could not write signal.\n", __FILE__, __LINE__);

  if ( fits_close_file(fptr, &status) )       /* close the FITS file */
    fprintf(stderr, "%s (%d): Could not close file.\n",
        __FILE__, __LINE__);

}


void s2let_fits_mw_read_spin_maps(double* fQ, double* fU, char* filename, int L)
{
  long     naxes, *naxis, npix, npercol, irow, Lread;
  int      status, hdutype, nfound, anynul;
  float    nulval;
  char     comment[FLEN_COMMENT];
  fitsfile *fptr;
  status = 0;

  if ( fits_open_file(&fptr, filename, READONLY, &status) )
    printerror( status );

  if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) )
    printerror( status );

  if (hdutype != BINARY_TBL)
    fprintf(stderr, "%s (%d): Extension is not binary!\n", __FILE__, __LINE__);

  if ( fits_read_key_lng(fptr, "NAXIS", &naxes, comment, &status) )
    printerror( status );

  naxis = (long *)malloc(((size_t)naxes)*sizeof(long));
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, naxes, naxis, &nfound, &status)
       || nfound != naxes )
    printerror( status );

  if ( fits_read_key_lng(fptr, "L", &Lread, comment, &status) )
    printerror(status);

  if( Lread != L )
    printf("Attention : read L = %li but you specified L = %i\n",Lread,L);

  npix = L * (2*L-1);
  printf("naxis[0] = %li\n",naxis[0]);
  printf("naxis[1] = %li\n",naxis[1]);

  if ( (npix%naxis[1]) != 0 )
    fprintf(stderr, "%s (%d): Problem with npix.\n", __FILE__, __LINE__);

  npercol = npix/naxis[1];
  nulval = -1.6375e30;
  for (irow = 0; irow < naxis[1]; irow++) {
    if ( fits_read_col(fptr, TDOUBLE, 1, irow+1, 1, npercol, &nulval,
           &(fQ[irow*npercol]), &anynul, &status) ) {
      printerror(status);
    }
  }
  for (irow = 0; irow < naxis[1]; irow++) {
    if ( fits_read_col(fptr, TDOUBLE, 2, irow+1, 1, npercol, &nulval,
           &(fU[irow*npercol]), &anynul, &status) ) {
      printerror(status);
    }
  }
  if ( fits_close_file(fptr, &status) )
    printerror( status );

}
void s2let_fits_mwss_read_spin_maps(double* fQ, double* fU, char* filename, int L)
{
  long     naxes, *naxis, npix, npercol, irow, Lread;
  int      status, hdutype, nfound, anynul;
  float    nulval;
  char     comment[FLEN_COMMENT];
  fitsfile *fptr;
  status = 0;

  if ( fits_open_file(&fptr, filename, READONLY, &status) )
    printerror( status );

  if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) )
    printerror( status );

  if (hdutype != BINARY_TBL)
    fprintf(stderr, "%s (%d): Extension is not binary!\n", __FILE__, __LINE__);

  if ( fits_read_key_lng(fptr, "NAXIS", &naxes, comment, &status) )
    printerror( status );

  naxis = (long *)malloc(((size_t)naxes)*sizeof(long));
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, naxes, naxis, &nfound, &status)
       || nfound != naxes )
    printerror( status );

  if ( fits_read_key_lng(fptr, "L", &Lread, comment, &status) )
    printerror(status);

  if( Lread != L )
    printf("Attention : read L = %li but you specified L = %i\n",Lread,L);

  long mwss_read;
  if ( fits_read_key_lng(fptr, "MWSS", &mwss_read, comment, &status) )
    printerror(status);

  if( mwss_read != 1 )
    printf("Attention : the file is not stored with MWSS sampling.\n");

  npix = (L+1) * 2*L;
  printf("naxis[0] = %li\n",naxis[0]);
  printf("naxis[1] = %li\n",naxis[1]);

  if ( (npix%naxis[1]) != 0 )
    fprintf(stderr, "%s (%d): Problem with npix.\n", __FILE__, __LINE__);

  npercol = npix/naxis[1];
  nulval = -1.6375e30;
  for (irow = 0; irow < naxis[1]; irow++) {
    if ( fits_read_col(fptr, TDOUBLE, 1, irow+1, 1, npercol, &nulval,
           &(fQ[irow*npercol]), &anynul, &status) ) {
      printerror(status);
    }
  }
  for (irow = 0; irow < naxis[1]; irow++) {
    if ( fits_read_col(fptr, TDOUBLE, 2, irow+1, 1, npercol, &nulval,
           &(fU[irow*npercol]), &anynul, &status) ) {
      printerror(status);
    }
  }
  if ( fits_close_file(fptr, &status) )
    printerror( status );

}

/*
void s2let_fits_mw_write_wav_maps(char* filename, complex double *f_wav,
    complex double *f_scal, int nmaps,
    int L, int N, int B, int J_min, int J, int multires) // TODO: ADD PARAMETERS STRUCT
{
  fitsfile *fptr;
  int i, status, hdutype;
  long firstrow, firstelem;
  int bitpix   =  SHORT_IMG;
  long naxis   =   0;
  long naxes[] = {0,0};
  int tfields   = 1;
  char extname[] = "BINTABLE";
  char *ttype[] = { "SIGNAL" };
  char *tform[] = { "1D" }; COMPLEX
  char *tunit[] = { " " };

  long npix = L * (2*L-1);
  status = 0;

  if (fits_create_file(&fptr, filename, &status))
    fprintf(stderr, "%s (%d): Could not create new fits file.\n",
      __FILE__, __LINE__);

  if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
    fprintf(stderr, "%s (%d): Could not create new image file.\n",
      __FILE__, __LINE__);

  if ( fits_movabs_hdu(fptr, 1, &hdutype, &status) )
    fprintf(stderr, "%s (%d): Could not move to first HDU.\n",
      __FILE__, __LINE__);


  for(i=0; i<nmaps; ++i){

    int bl = L;

    if ( fits_create_tbl( fptr, BINARY_TBL, npix, tfields, ttype, tform,
        tunit, extname, &status) )
      fprintf(stderr, "%s (%d): Could not create new binary table.\n",
        __FILE__, __LINE__);

    if (fits_write_key(fptr, TINT, "L", &bl,
           "Resolution parameter", &status))
      fprintf(stderr, "%s (%d): Could not write L keyword.\n",
        __FILE__, __LINE__);

    firstrow  = 1;
    firstelem = 1;
    npix = bl * (2*bl-1);
    if (fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, npix, fwav[, &status))
      fprintf(stderr, "%s (%d): Could not write signal.\n", __FILE__, __LINE__);

  }

  if ( fits_close_file(fptr, &status) )
    fprintf(stderr, "%s (%d): Could not close file.\n",
      __FILE__, __LINE__);

}


void s2let_fits_mw_read_wav_maps(complex double *f_wav,
    complex double *f_scal, char* filename, int nmaps,
    int L, int N, int B, int J_min, int J, int multires) // TODO: ADD PARAMETERS STRUCT
{
  long     naxes, *naxis, npix, npercol, irow, Lread;
  int      status, hdutype, nfound, anynul;
  float    nulval;
  char     comment[FLEN_COMMENT];
  fitsfile *fptr;
  status = 0;

  if ( fits_open_file(&fptr, filename, READONLY, &status) )
    printerror( status );

  if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) )
    printerror( status );

  if (hdutype != BINARY_TBL)
    fprintf(stderr, "%s (%d): Extension is not binary!\n", __FILE__, __LINE__);

  if ( fits_read_key_lng(fptr, "NAXIS", &naxes, comment, &status) )
    printerror( status );

  naxis = (long *)malloc(((size_t)naxes)*sizeof(long));
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, naxes, naxis, &nfound, &status)
       || nfound != naxes )
    printerror( status );

  if ( fits_read_key_lng(fptr, "L", &Lread, comment, &status) )
    printerror(status);

  if( Lread != L )
    printf("Attention : read L = %li but you specified L = %i\n",Lread,L);

  npix = L * (2*L-1);
  printf("naxis[0] = %i\n",naxis[0]);
  printf("naxis[1] = %i\n",naxis[1]);

  if ( (npix%naxis[1]) != 0 )
    fprintf(stderr, "%s (%d): Problem with npix.\n", __FILE__, __LINE__);

  npercol = npix/naxis[1];
  nulval = -1.6375e30;
  for (irow = 0; irow < naxis[1]; irow++) {
    if ( fits_read_col(fptr, TDOUBLE, 1, irow+1, 1, npercol, &nulval,
           &(fQ[irow*npercol]), &anynul, &status) ) {
      printerror(status);
    }
  }
  for (irow = 0; irow < naxis[1]; irow++) {
    if ( fits_read_col(fptr, TDOUBLE, 2, irow+1, 1, npercol, &nulval,
           &(fU[irow*npercol]), &anynul, &status) ) {
      printerror(status);
    }
  }
  if ( fits_close_file(fptr, &status) )
    printerror( status );

}*/
