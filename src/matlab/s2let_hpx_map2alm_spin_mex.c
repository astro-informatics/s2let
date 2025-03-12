// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include <s2let.h>
#include "mex.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))

/**
 * MATLAB interface: s2let_hpx_map2alm_spin_mex.
 * This function for internal use only.
 *
 * Usage:
 *   alm!, almU = s2let_hpx_map2alm_spin_mex(fQ, fU, nside, L, spin);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int i, nside, L, spin, f_m, f_n;
  complex double *flmQ = NULL, *flmU = NULL;
  double *f_real = NULL, *flm_real = NULL, *flm_imag = NULL, *fQ_r = NULL, *fU_r = NULL;
  int iin = 0, iout = 0;


  // Parse input dataset f
  iin = 0;
  f_m = mxGetM(prhs[iin]);
  f_n = mxGetN(prhs[iin]);
  f_real = mxGetPr(prhs[iin]);
  fQ_r = (double*)malloc(f_m * f_n * sizeof(double));
  for (i=0; i<f_m*f_n; i++)
    fQ_r[i] = f_real[i];

  iin = 1;
  f_m = mxGetM(prhs[iin]);
  f_n = mxGetN(prhs[iin]);
  f_real = mxGetPr(prhs[iin]);
  fU_r = (double*)malloc(f_m * f_n * sizeof(double));
  for (i=0; i<f_m*f_n; i++)
    fU_r[i] = f_real[i];

  // Parse HEALPIX parameter nside
  iin = 2;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_hpx_map2alm_spin_mex:InvalidInput:healpixParameter",
          "HEALPIX parameter nside must be integer.");
  }
  nside = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)nside || nside <= 1)
    mexErrMsgIdAndTxt("s2let_hpx_map2alm_spin_mex:InvalidInput:healpixParameter",
          "Healpix parameter nside must be positive integer greater than 2");

  if( f_m*f_n != 12*nside*nside )
    mexErrMsgIdAndTxt("s2let_hpx_map2alm_spin_mex:InvalidInput:LbandLimit",
          "nside must correspond to the sampling scheme, i.e. f = 12*nside*nside samples.");

  // Parse harmonic band-limit L
  iin = 3;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_hpx_map2alm_spin_mex:InvalidInput:LbandLimit",
          "Harmonic band-limit L must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("s2let_hpx_map2alm_spin_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit L must be positive integer.");

  // Parse spin
  iin = 4;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_hpx_alm2map_spin_mex:InvalidInput:LbandLimit",
          "spin must be integer.");
  }
  spin = (int)mxGetScalar(prhs[iin]);

  // Perform harmonic transform
  s2let_allocate_lm(&flmQ, L);
  s2let_allocate_lm(&flmU, L);
  s2let_hpx_map2alm_spin_real(flmQ, flmU, fQ_r, fU_r, nside, L, spin);

  // Output flm's
  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(L * L, 1, mxCOMPLEX);
  flm_real = mxGetPr(plhs[iout]);
  flm_imag = mxGetPi(plhs[iout]);
  for (i=0; i<L*L; i++){
    flm_real[i] = creal( flmQ[i] );
    flm_imag[i] = cimag( flmQ[i] );
  }

  iout = 1;
  plhs[iout] = mxCreateDoubleMatrix(L * L, 1, mxCOMPLEX);
  flm_real = mxGetPr(plhs[iout]);
  flm_imag = mxGetPi(plhs[iout]);
  for (i=0; i<L*L; i++){
    flm_real[i] = creal( flmU[i] );
    flm_imag[i] = cimag( flmU[i] );
  }

  free(flmQ);
  free(fQ_r);
  free(flmU);
  free(fU_r);


}
