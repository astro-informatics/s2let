// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "s2let/s2let.h"

/*!
 * PROGRAM : s2let_transform_axisym_mw_analysis_real
 * COMMAND : bin/s2let_transform_axisym_mw_analysis_real file B J_min multires
 * ARGUMENTS :
 * - file : input MW map
 * - B : wavelet parameter
 * - J_min : first wavelet scale to use
 * - multires : multiresolution flag (1: activated, 0: off)
 * OUTPUT : fits files containing the wavelet MW maps
 */
int main(int argc, char *argv[]) {
  printf("--------------------------------------------------\n");
  printf("S2LET library : axisymmetric wavelet transform\n");
  printf("Real signal, MW sampling\n");
  printf("--------------------------------------------------\n");

  char file[100];
  if (sscanf(argv[1], "%s", file) != 1)
    exit(-2);
  printf("Input MW map : %s\n", file);
  const int L = s2let_fits_mw_read_bandlimit(file);
  printf("- Detected bandlimit L = %i\n", L);
  int B, J_min, multires;
  if (sscanf(argv[2], "%i", &B) != 1)
    exit(-2);
  if (sscanf(argv[3], "%i", &J_min) != 1)
    exit(-2);
  if (sscanf(argv[4], "%i", &multires) != 1)
    exit(-2);

  s2let_parameters_t parameters = {};

  parameters.B = B;
  parameters.L = L;
  parameters.J_min = J_min;

  printf("Parameters for wavelet transform :\n");
  int J = s2let_j_max(&parameters);
  printf("- Multiresolution flag : %i\n", multires);
  printf("- Wavelet parameter : %i\n", B);
  printf("- Total number of wavelets : %i\n", J);
  printf("- First wavelet scale to be used : %i\n", J_min);

  // Read MW map from file
  double *f = (double *)calloc(L * (2 * L - 1), sizeof(double));
  s2let_fits_mw_read_map(f, file, L);
  printf("File successfully read from file\n");

  printf("Performing wavelet decomposition...");
  fflush(NULL);
  double *f_wav, *f_scal;
  if (multires) {
    s2let_transform_axisym_allocate_mw_f_wav_multires_real(
        &f_wav, &f_scal, &parameters);
    s2let_transform_axisym_wav_analysis_mw_multires_real(f_wav, f_scal, f, &parameters);
  } else {
    s2let_transform_axisym_allocate_mw_f_wav_real(&f_wav, &f_scal, &parameters);
    s2let_transform_axisym_wav_analysis_mw_real(f_wav, f_scal, f, &parameters);
  }
  printf("done\n");

  // Output the wavelets to FITS files
  printf("Writing wavelet maps to FITS files\n");
  char outfile[100];
  char params[100];
  sprintf(params, "%d%s%d%s%d", L, "_", B, "_", J_min);
  int j, bl;      // Explicitly compute the maximum wavelet scale
  int offset = 0; // Start with the first wavelet
  char fileroot[100];
  sscanf(file, "%[^.]", fileroot);
  printf("File root = %s\n", fileroot);
  for (j = J_min; j <= J; j++) {
    sprintf(outfile, "%s%s%s%s%d%s", fileroot, "_wav_", params, "_", j, ".fits");
    printf("- Outfile_wav[j=%i] = %s\n", j, outfile);
    remove(outfile); // In case the file exists
    if (multires)
      bl = MIN(s2let_bandlimit(j, &parameters), L);
    else
      bl = L;
    s2let_fits_mw_write_map(
        outfile, f_wav + offset, bl); // Now write the map to fits file
    offset += (2 * bl - 1) * bl;      // Go to the next wavelet
  }
  // Finally write the scaling function
  sprintf(outfile, "%s%s%s%s", fileroot, "_scal_", params, ".fits");
  printf("- Outfile_scal = %s\n", outfile);
  remove(outfile); // In case the file exists
  if (multires)
    bl = MIN(s2let_bandlimit(J_min - 1, &parameters), L);
  else
    bl = L;
  s2let_fits_mw_write_map(outfile, f_scal, bl); // Now write the map to fits file

  printf("--------------------------------------------------\n");

  return 0;
}
