// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include <s2let.h>
#include "mex.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))

/**
 * MATLAB interface: s2let_transform_axisym_analysis_hpx.
 * This function for internal use only.
 * Compute axisymmetric wavelet transform (analysis)
 * with output in pixel space. Format : HEALPIX maps.
 *
 * Usage:
 *   [f_wav, f_scal] = ...
 *        s2let_transform_axisym_analysis_hpx_mex(f, nside, B, L, J_min);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int i, nside, L, J_min, f_m, f_n;
  double B;
  s2let_parameters_t parameters = {};
  double *f_wav_real, *f_scal_real, *f_real;
  double *f_wav_r = NULL, *f_scal_r = NULL, *f_r = NULL;
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=5) {
    mexErrMsgIdAndTxt("s2let_transform_axisym_analysis_hpx_mex:InvalidInput:nrhs",
          "Require six inputs.");
  }
  if(nlhs!=2) {
    mexErrMsgIdAndTxt("s2let_transform_axisym_analysis_hpx_mex:InvalidOutput:nlhs",
          "Require two outputs.");
  }
  // Parse input dataset f
  iin = 0;
  f_m = mxGetM(prhs[iin]);
  f_n = mxGetN(prhs[iin]);
  f_real = mxGetPr(prhs[iin]);
  f_r = (double*)malloc(f_m * f_n * sizeof(double));
  for (i=0; i<f_m*f_n; i++)
    f_r[i] = f_real[i];

  // Parse HEALPIX parameter nside
  iin = 1;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_axisym_analysis_hpx_mex:InvalidInput:healpixParameter",
          "HEALPIX parameter nside must be integer.");
  }
  nside = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)nside || nside <= 1)
    mexErrMsgIdAndTxt("s2let_transform_axisym_analysis_hpx_mex:InvalidInput:healpixParameter",
          "Healpix parameter nside must be positive integer greater than 2");

  if( f_m*f_n != 12*nside*nside )
    mexErrMsgIdAndTxt("s2let_transform_axisym_analysis_hpx_mex:InvalidInput:LbandLimit",
          "nside must correspond to the sampling scheme, i.e. f = 12*nside*nside samples.");

  // Parse wavelet parameter B
  iin = 2;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_axisym_analysis_hpx_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B must be integer.");
  }
  B = (double)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)B || B <= 1)
    mexErrMsgIdAndTxt("s2let_transform_axisym_analysis_hpx_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B must be positive integer greater than 2");

  // Parse harmonic band-limit L
  iin = 3;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_axisym_analysis_hpx_mex:InvalidInput:LbandLimit",
          "Harmonic band-limit L must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("s2let_transform_axisym_analysis_hpx_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit L must be positive integer.");

  // Parse first scale J_min
  iin = 4;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_axisym_analysis_hpx_mex:InvalidInput:Jmin",
          "First scale J_min must be integer.");
  }
  J_min = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)J_min || J_min < 0)
    mexErrMsgIdAndTxt("s2let_transform_axisym_analysis_hpx_mex:InvalidInput:Jmin",
          "First scale J_min must be positive integer.");

  parameters.B = B;
  parameters.L = L;
  parameters.J_min = J_min;

  // Compute ultimate scale J_max
  int J = s2let_j_max(&parameters);

  if( J_min > J+1 ) {
    mexErrMsgIdAndTxt("s2let_transform_axisym_analysis_hpx_mex:InvalidInput:Jmin",
          "First scale J_min must be larger than that!");
  }

  // Perform wavelet transform in harmonic space and then reconstruction.
  s2let_transform_axisym_allocate_hpx_f_wav_real(&f_wav_r, &f_scal_r, nside, &parameters);
  s2let_transform_axisym_wav_analysis_hpx_real(f_wav_r, f_scal_r, f_r, nside, &parameters);

  // Compute size of wavelet array
  int wavsize, scalsize;
  wavsize = (J+1-J_min) * 12 * nside * nside;
  scalsize = 12 * nside * nside;

  // Output wavelets
  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(1, wavsize, mxREAL);
  f_wav_real = mxGetPr(plhs[iout]);
  for (i=0; i<wavsize; i++){
    f_wav_real[i] = creal(f_wav_r[i]);
  }

  iout = 1;
  plhs[iout] = mxCreateDoubleMatrix(1, scalsize, mxREAL);
  f_scal_real = mxGetPr(plhs[iout]);
  for (i=0; i<scalsize; i++)
    f_scal_real[i] = creal(f_scal_r[i]);


  free(f_r);
  free(f_wav_r);
  free(f_scal_r);


}
