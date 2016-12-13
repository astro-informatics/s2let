// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include <s2let.h>
#include <s2let_mex.h>
#include <string.h>
#include "mex.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))

/**
 * MATLAB interface: s2let_transform_analysis_lm2wav_mex.
 * This function for internal use only.
 * Compute spin directional wavelet transform (analysis)
 * with input in harmonic space and output in pixel space.
 *
 * Usage:
 *   [f_wav, f_scal] = ...
 *        s2let_transform_analysis_lm2wav_mex(flm, B, L, J_min, N, spin, reality, upsample,
 *                                            spin_lowered, original_spin, sampling_scheme);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  int i, j, B, L, J_min, N, spin, flm_m, flm_n, flm_size, reality, upsample, normalization, original_spin;
  char sampling_str[S2LET_STRING_LEN];
  s2let_sampling_t sampling_scheme;
  s2let_parameters_t parameters = {};
  double *f_wav_real, *f_scal_real, *flm_real;
  double *f_wav_imag, *f_scal_imag, *flm_imag;
  complex double *f_wav = NULL, *f_scal = NULL, *flm= NULL;
  double *f_wav_r = NULL, *f_scal_r = NULL;
  int iin = 0, iout = 0;
  // Check number of arguments
  if(nrhs!=11) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_lm2wav_mex:InvalidInput:nrhs",
          "Require eleven inputs.");
  }
  if(nlhs!=2) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_lm2wav_mex:InvalidOutput:nlhs",
          "Require two outputs.");
  }

  // Parse reality flag
  iin = 6;
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("s2let_transform_analysis_lm2wav_mex:InvalidInput:reality",
          "Reality flag must be logical.");
  reality = mxIsLogicalScalarTrue(prhs[iin]);

  // Parse multiresolution flag
  iin = 7;
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("s2let_transform_analysis_lm2wav_mex:InvalidInput:upsample",
          "Multiresolution flag must be logical.");
  upsample = mxIsLogicalScalarTrue(prhs[iin]);

  /* Parse sampling scheme method. */
  iin = 10;
  if( !mxIsChar(prhs[iin]) ) {
      mexErrMsgIdAndTxt("s2let_transform_analysis_lm2wav_mex:InvalidInput:samplingSchemeChar",
                        "Sampling scheme must be string.");
  }
  int len = (mxGetM(prhs[iin]) * mxGetN(prhs[iin])) + 1;
  if (len >= S2LET_STRING_LEN)
      mexErrMsgIdAndTxt("s2let_transform_analysis_lm2wav_mex:InvalidInput:samplingSchemeTooLong",
                        "Sampling scheme exceeds string length.");
  mxGetString(prhs[iin], sampling_str, len);

  if (strcmp(sampling_str, S2LET_SAMPLING_MW_STR) == 0)
      sampling_scheme = S2LET_SAMPLING_MW;
  else if (strcmp(sampling_str, S2LET_SAMPLING_MW_SS_STR) == 0)
      sampling_scheme = S2LET_SAMPLING_MW_SS;
  else
      mexErrMsgIdAndTxt("s2let_transform_analysis_lm2wav_mex:InvalidInput:samplingScheme",
                        "Invalid sampling scheme.");

  // Parse normalization flag
  iin = 8;
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("s2let_transform_analysis_lm2wav_mex:InvalidInput:spinlowered",
          "SpinLowered flag must be logical.");
  if (mxIsLogicalScalarTrue(prhs[iin]))
    normalization = S2LET_WAV_NORM_SPIN_LOWERED;
  else
    normalization = S2LET_WAV_NORM_DEFAULT;

  // Parse original spin
  iin = 9;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_lm2wav_mex:InvalidInput:spinloweredfrom",
          "SpinLoweredFrom must be integer.");
  }
  original_spin = (int)mxGetScalar(prhs[iin]);

  // Parse input harmonic coefficients flm
  iin = 0;
  flm_m = mxGetM(prhs[iin]);
  flm_n = mxGetN(prhs[iin]);
  flm_size = flm_m * flm_n;
  flm = malloc(flm_size * sizeof(*flm));
  flm_real = mxGetPr(prhs[iin]);
  flm_imag = mxGetPi(prhs[iin]);

  for (i = 0; i < flm_size; ++i)
    flm[i] = flm_real[i] + I * flm_imag[i];

  // Parse wavelet parameter B
  iin = 1;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_lm2wav_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B must be integer.");
  }
  B = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)B || B <= 1)
    mexErrMsgIdAndTxt("s2let_transform_analysis_lm2wav_mex:InvalidInput:waveletParameter",
          "Wavelet parameter B must be positive integer greater than 2");

  // Parse harmonic band-limit L
  iin = 2;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_lm2wav_mex:InvalidInput:LbandLimit",
          "Harmonic band-limit L must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);

  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("s2let_transform_analysis_lm2wav_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit L must be positive integer.");

  if( flm_size != L*L ) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_lm2wav_mex:InvalidInput:LbandLimit",
          "L must correspond to input size, i.e. L*L harmonic coefficients.");
  }

  // Parse first scale J_min
  iin = 3;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_lm2wav_mex:InvalidInput:Jmin",
          "First scale J_min must be integer.");
  }
  J_min = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)J_min || J_min < 0)
    mexErrMsgIdAndTxt("s2let_transform_analysis_lm2wav_mex:InvalidInput:Jmin",
          "First scale J_min must be positive integer.");

  parameters.B = B;
  parameters.L = L;
  parameters.J_min = J_min;

  // Compute ultimate scale J_max
  int J = s2let_j_max(&parameters);

  if( J_min > J+1 ) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_lm2wav_mex:InvalidInput:Jmin",
          "First scale J_min must be larger than that!");
  }

    // Parse azimuthal/directional band-limit N
  iin = 4;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_lm2wav_mex:InvalidInput:NbandLimit",
          "Azimuthal/directional band-limit N must be integer.");
  }
  N = (int)mxGetScalar(prhs[iin]);

  // Parse spin
  iin = 5;
  if( !mxIsDouble(prhs[iin]) ||
      mxIsComplex(prhs[iin]) ||
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("s2let_transform_analysis_lm2wav_mex:InvalidInput:spin",
          "spin must be integer.");
  }
  spin = (int)mxGetScalar(prhs[iin]);

  if (reality && spin)
    mexErrMsgIdAndTxt("s2let_transform_analysis_lm2wav_mex:InvalidInput:realspin",
                      "Real signals must have spin zero.");

  parameters.N = N;
  parameters.spin = spin;
  parameters.upsample = upsample;
  parameters.normalization = normalization;
  parameters.original_spin = original_spin;
  parameters.reality = reality;
  parameters.sampling_scheme = sampling_scheme;

  // Perform wavelet transform in harmonic space and then reconstruction.
  if(reality){
      s2let_allocate_f_wav_real(&f_wav_r, &f_scal_r, &parameters);
      s2let_analysis_lm2wav_real(f_wav_r, f_scal_r, flm, &parameters);
  }else{
      s2let_allocate_f_wav(&f_wav, &f_scal, &parameters);
      s2let_analysis_lm2wav(f_wav, f_scal, flm, &parameters);
  }

  // Compute size of wavelet array
  int bandlimit, wavsize = 0, scalsize = 0;
  so3_parameters_t so3_parameters = {};
  so3_parameters.N = N;
  so3_parameters.sampling_scheme = sampling_scheme;
  so3_parameters.steerable = 0; // When direct routines support steerable so3 change this to 1.
  if(upsample){
    so3_parameters.L = L;
    wavsize = (J+1-J_min) * so3_sampling_f_size(&so3_parameters);
    if (sampling_scheme == S2LET_SAMPLING_MW_SS)
        scalsize = (L+1) * 2*L;
    else
        scalsize = L * (2*L - 1);
  }else{
    for (j = J_min; j <= J; j++){
        bandlimit = MIN(s2let_bandlimit(j, &parameters), L);
        so3_parameters.L = bandlimit;
        wavsize += so3_sampling_f_size(&so3_parameters);
     }
     bandlimit = MIN(s2let_bandlimit(J_min-1, &parameters), L);
    if (sampling_scheme == S2LET_SAMPLING_MW_SS)
        scalsize = (bandlimit+1) * 2 * bandlimit;
    else
        scalsize = bandlimit * (2 * bandlimit - 1);
  }

  // Output wavelets
  if (reality)
  {
    iout = 0;
    plhs[iout] = mxCreateDoubleMatrix(1, wavsize, mxREAL);
    f_wav_real = mxGetPr(plhs[iout]);
    for (i=0; i<wavsize; i++){
      f_wav_real[ i ] = f_wav_r[ i ];
    }

    iout = 1;
    plhs[iout] = mxCreateDoubleMatrix(1, scalsize, mxREAL);
    f_scal_real = mxGetPr(plhs[iout]);
    for (i=0; i<scalsize; i++){
      f_scal_real[i] = f_scal_r[i];
    }
  }
  else
  {
    iout = 0;
    plhs[iout] = mxCreateDoubleMatrix(1, wavsize, mxCOMPLEX);
    f_wav_real = mxGetPr(plhs[iout]);
    f_wav_imag = mxGetPi(plhs[iout]);
    for (i=0; i<wavsize; i++){
      f_wav_real[ i ] = creal( f_wav[ i ] );
      f_wav_imag[ i ] = cimag( f_wav[ i ] );
    }

    iout = 1;
    plhs[iout] = mxCreateDoubleMatrix(1, scalsize, mxCOMPLEX);
    f_scal_real = mxGetPr(plhs[iout]);
    f_scal_imag = mxGetPi(plhs[iout]);
    for (i=0; i<scalsize; i++){
      f_scal_real[i] = creal( f_scal[i] );
      f_scal_imag[i] = cimag( f_scal[i] );
    }
  }

  free(flm);
  if(reality){
    free(f_wav_r);
    free(f_scal_r);
  }else{
    free(f_wav);
    free(f_scal);
  }

}
