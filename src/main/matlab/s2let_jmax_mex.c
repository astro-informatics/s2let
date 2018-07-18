// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include <s2let.h>
#include "mex.h"

/**
 * MATLAB interface: s2let_jmax.
 * This function for internal use only.
 * Compute maximum wavelet scale.
 *
 * Usage:
 *   bl = s2let_jmax(L, B);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  int L;
  double B;
  s2let_parameters_t parameters = {};
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=2) {
    mexErrMsgIdAndTxt("s2let_axisym_tiling_mex:InvalidInput:nrhs",
		      "Require three inputs.");
  }
  if(nlhs!=1) {
    mexErrMsgIdAndTxt("s2let_axisym_tiling_mex:InvalidOutput:nlhs",
		      "Require one outputs.");
  }

  // Parse harmonic band-limit L
  iin = 0;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_axisym_tiling_mex:InvalidInput:LbandLimit",
		      "Harmonic band-limit L must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);

  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("s2let_axisym_tiling_mex:InvalidInput:bandLimitNonInt",
		      "Harmonic band-limit L must be positive integer.");

  // Parse wavelet parameter B
  iin = 1;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_axisym_tiling_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B must be integer.");
  }
  B = (double)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)B || B <= 1)
    mexErrMsgIdAndTxt("s2let_axisym_tiling_mex:InvalidInput:bandLimitNonInt",
          "Wavelet parameter B must be positive integer greater than 2");

  if( B >= L ) {
    mexErrMsgIdAndTxt("s2let_axisym_tiling_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B must be smaller than L!");
  }

  parameters.B = B;
  parameters.L = L;

  // Compute ultimate scale J_max
  int Jmax = s2let_j_max(&parameters);
  double *out;

  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(1, 1, mxREAL);
  out = mxGetPr(plhs[iout]);
  out[0] = Jmax;

}
