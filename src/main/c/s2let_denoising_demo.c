// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <ssht/ssht.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "s2let/s2let.h"

void s2let_lm_random_flm_real_sigma(
    complex double *flm, int L, int seed, double sigmanoise) {
  int el, m, msign, i, i_op;
  for (el = 0; el < L; el++) {
    m = 0;
    i = el * el + el + m;
    flm[i] = (2.0 * ran2_dp(seed) - 1.0);
    for (m = 1; m <= el; m++) {
      i = el * el + el + m;
      flm[i] = sigmanoise * (2.0 * ran2_dp(seed) - 1.0) +
               I * sigmanoise * (2.0 * ran2_dp(seed) - 1.0);
      i_op = el * el + el - m;
      msign = m & 1;
      msign = 1 - msign - msign; // (-1)^m
      flm[i_op] = msign * conj(flm[i]);
    }
  }
}

double waveletpower(complex double *wav_lm, int L) {
  int i;
  double totalpower = 0;
  for (i = 0; i < L * L; i++)
    totalpower += wav_lm[i] * conj(wav_lm[i]);
  return totalpower;
}

void hard_threshold_real(
    double *g_wav, const double *threshold, const s2let_parameters_t *parameters) {
  int L = parameters->L;
  int J_min = parameters->J_min;
  int N = parameters->N;

  int J = s2let_j_max(parameters);
  int i, j, offset = 0;
  for (j = J_min; j <= J; j++) {
    int bl = parameters->upsample ? L : MIN(s2let_bandlimit(j, parameters), L);
    for (i = 0; i < N * bl * (2 * bl - 1); i++) {
      if (cabs(g_wav[offset + i]) < threshold[j - J_min])
        g_wav[offset + i] = 0;
    }
    offset += N * bl * (2 * bl - 1);
  }
}

/*!
 * PROGRAM : s2let_denoising_demo
 * COMMAND : bin/s2let_denoising_demo
 * ARGUMENTS : none
 */
int main(int argc, char *argv[]) {
  s2let_parameters_t parameters = {};

  const int seed = (int)(10000.0 * (double)clock() / (double)CLOCKS_PER_SEC);
  // PARAMETERS
  const double SNR_in = 10.0; // Input SNR
  const int nsigma = 3;       // Number of sigmas for hard thresholding
  const int upsample = 0;     // Multiresolution flag
  const double B = 2;         // Wavelet parameters
  const int N = 4;            // Azimuthal band-limit
  const int J_min = 0;        // First wavelet scale to use

  char outfile[100];
  int i, j;
  double *f, *noise, *g, *g_wav, *g_scal, *scal_l, *f_denoised, *remaining_noise;
  complex double *noise_lm, *wav_lm;

  parameters.B = B;
  parameters.J_min = J_min;
  parameters.N = N;
  parameters.upsample = upsample;
  parameters.reality = 1;

  printf("--------------------------------------------------\n");
  printf(" S2LET library : denoising example\n");
  printf(" Earth tomography signal, MW sampling\n");
  printf("--------------------------------------------------\n");

  char file[100] = "data/earth_tomo_mw_128.fits";
  printf(" Reading file %s\n", file);
  const int L = s2let_fits_mw_read_bandlimit(file);
  parameters.L = L;
  printf(" - Detected bandlimit L = %i\n", L);
  int J = s2let_j_max(&parameters);
  printf(" Parameters for wavelet denoising :\n");
  s2let_switch_wavtype(1);
  printf(" - Input SNR : %f\n", SNR_in);
  printf(" - Sigma threshold : %i\n", nsigma);
  printf(" - upsample flag : %i\n", upsample);
  printf(" - Wavelet parameter : %i\n", B);
  printf(" - Total number of wavelets : %i\n", J);
  printf(" - First wavelet scale to be used : %i\n", J_min);

  s2let_allocate_mw_real(&f, L);
  s2let_fits_mw_read_map(f, file, L); // Read MW map from file
  for (i = 0; i < L * (2 * L - 1); i++)
    f[i] = f[i] * 1000;
  printf(" File successfully read from file\n");

  // Compute noise standard deviation and generate noise
  double sigmanoise = sqrt(pow(10.0, -SNR_in / 10.0) * s2let_mw_power_real(f, L));
  printf(" - Std dev of the noise (to match SNR) = %f\n", sigmanoise);
  s2let_allocate_lm(&noise_lm, L);
  s2let_lm_random_flm_real_sigma(noise_lm, L, seed, sigmanoise);
  double SNR_actual =
      10.0 * log10(s2let_mw_power_real(f, L) / s2let_lm_power(noise_lm, L));
  printf(" - Actual (realised) SNR = %f\n", SNR_actual);

  // Add noise to the signal in real space
  printf(" Contaminating the signal with this noise...");
  fflush(NULL);
  s2let_allocate_mw_real(&noise, L);
  s2let_mw_alm2map_real(noise, noise_lm, L);
  s2let_allocate_mw_real(&g, L);
  for (i = 0; i < L * (2 * L - 1); i++)
    g[i] = f[i] + noise[i];
  printf(" done\n");

  printf(" Performing wavelet decomposition...");
  fflush(NULL);
  // Perform wavelet analysis from scratch with all signals given as MW maps
  s2let_allocate_f_wav_real(&g_wav, &g_scal, &parameters);
  s2let_analysis_px2wav_real(g_wav, g_scal, g, &parameters);
  printf(" done\n");

  // Compute simple threshold for needlet coefficients based on noise model
  printf(" Construct the threshold rule for the Gaussian noise\n");
  s2let_tiling_wavelet_allocate(&wav_lm, &scal_l, &parameters);
  s2let_tiling_wavelet(wav_lm, scal_l, &parameters);
  double *threshold = (double *)calloc((J - J_min + 1), sizeof(double));
  for (j = J_min; j <= J; j++)
    threshold[j - J_min] =
        sigmanoise * nsigma * sqrt(waveletpower(wav_lm + j * L * L, L));

  printf(" Hard thresholding the wavelets...");
  fflush(NULL);
  s2let_allocate_mw_real(&f_denoised, L);
  hard_threshold_real(g_wav, threshold, &parameters);
  s2let_synthesis_wav2px_real(f_denoised, g_wav, g_scal, &parameters);
  printf(" done\n");

  // Remaining noise
  s2let_allocate_mw_real(&remaining_noise, L);
  for (i = 0; i < L * (2 * L - 1); i++)
    remaining_noise[i] = f_denoised[i] - f[i];

  // SNR after denoising
  double SNR_denoised =
      10.0 * log10(s2let_mw_power_real(f, L) / s2let_mw_power_real(remaining_noise, L));
  printf(" -> SNR before denoising = %f\n", SNR_actual);
  printf(" -> SNR after denoising  = %f\n", SNR_denoised);

  // Finally write the denoised signal
  printf(" Write output files\n");
  sprintf(outfile, "%s%s%s", "data/real_signal", "_U_input", ".fits");
  printf(" Outfile = %s\n", outfile);
  remove(outfile);                        // In case the file exists
  s2let_fits_mw_write_map(outfile, f, L); // Now write the map to fits file
  sprintf(outfile, "%s%s%s", "data/real_signal", "_U_noise", ".fits");
  printf(" Outfile = %s\n", outfile);
  remove(outfile);                            // In case the file exists
  s2let_fits_mw_write_map(outfile, noise, L); // Now write the map to fits file
  sprintf(outfile, "%s%s%s", "data/real_signal", "_U_input_noise", ".fits");
  printf(" Outfile = %s\n", outfile);
  remove(outfile);                        // In case the file exists
  s2let_fits_mw_write_map(outfile, g, L); // Now write the map to fits file
  sprintf(outfile, "%s%s%s", "data/real_signal", "_U_denoised", ".fits");
  printf(" Outfile = %s\n", outfile);
  remove(outfile);                                 // In case the file exists
  s2let_fits_mw_write_map(outfile, f_denoised, L); // Now write the map to fits file

  free(f);
  free(noise);
  free(g);
  free(g_wav);
  free(g_scal);
  free(scal_l);
  free(f_denoised);
  free(remaining_noise);
  free(noise_lm);

  printf("--------------------------------------------------\n");
  return 0;
}
