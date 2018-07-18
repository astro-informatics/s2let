// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include <s2let.h>
#include "mex.h"

/**
 * MATLAB interface: s2let_axisym_tiling_mex.
 * This function for internal use only.
 * Compute tiling in harmonic space for axisymmetric wavelets.
 *
 * Usage:
 *   [kappa kappa0] = s2let_axisym_tiling_mex(B, L, J_min);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  int L, J_min;
  double B;
  s2let_parameters_t parameters = {};
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=3) {
    mexErrMsgIdAndTxt("s2let_axisym_tiling_mex:InvalidInput:nrhs",
		      "Require three inputs.");
  }
  if(nlhs!=2) {
    mexErrMsgIdAndTxt("s2let_axisym_tiling_mex:InvalidOutput:nlhs",
		      "Require two outputs.");
  }

  // Parse wavelet parameter B
  iin = 0;
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

  // Parse harmonic band-limit L
  iin = 1;
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

  if( B >= L ) {
    mexErrMsgIdAndTxt("s2let_axisym_tiling_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B must be smaller than L!");
  }

  // Parse first scale J_min
  iin = 2;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_axisym_tiling_mex:InvalidInput:Jmin",
          "First scale J_min must be integer.");
  }
  J_min = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)J_min || J_min < 0)
    mexErrMsgIdAndTxt("s2let_axisym_tiling_mex:InvalidInput:Jmin",
          "First scale J_min must be positive integer.");

  parameters.B = B;
  parameters.L = L;
  parameters.J_min = J_min;

  // Compute ultimate scale J_max
  int J = s2let_j_max(&parameters);

  if( J_min > J+1 ) {
    mexErrMsgIdAndTxt("s2let_axisym_tiling_mex:InvalidInput:Jmin",
          "First scale J_min must be larger than that!");
  }

  // Allocate arrays
  double *kappa = (double*)calloc((J+1) * L, sizeof(double));
  double *kappa0 = (double*)calloc(L, sizeof(double));

  // Run S2LET function
  s2let_tiling_axisym(kappa, kappa0, &parameters);


  // Output kappa and kappa0
  double *kappa_out, *kappa0_out;
  int l, j;

  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(J+1, L, mxREAL);
  kappa_out = mxGetPr(plhs[iout]);
  for(j=0; j<=J; j++)
    for(l=0; l<L; l++)
      kappa_out[l*(J+1)+j] = kappa[l+j*L];

  iout = 1;
  plhs[iout] = mxCreateDoubleMatrix(1, L, mxREAL);
  kappa0_out = mxGetPr(plhs[iout]);
  for(l=0; l<L; l++)
    kappa0_out[l] = kappa0[l];

  free(kappa);
  free(kappa0);

}
