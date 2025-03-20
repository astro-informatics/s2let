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
 * PROGRAM : s2let_transform_axisym_mw_synthesis_real
 * COMMAND : bin/s2let_transform_axisym_mw_synthesis_real file B J_min L
 * ARGUMENTS :
 * - file : root for input MW wavelet maps
 * - B : wavelet parameter
 * - J_min : first wavelet scale to use
 * - L : band-limit of the output map / of the finest wavelet map
 * OUTPUT : fits file containing the reconstructed MW map
 */
int main(int argc, char *argv[]) {
  printf("--------------------------------------------------\n");
  printf("S2LET library : axisymmetric wavelet transform\n");
  printf("Real signal, MW sampling\n");
  printf("--------------------------------------------------\n");

  char fileroot[100];
  if (sscanf(argv[1], "%s", fileroot) != 1)
    exit(-2);
  int multires, B, J_min, L;
  if (sscanf(argv[2], "%i", &B) != 1)
    exit(-2);
  if (sscanf(argv[3], "%i", &J_min) != 1)
    exit(-2);
  if (sscanf(argv[4], "%i", &L) != 1)
    exit(-2);

  s2let_parameters_t parameters = {};

  parameters.B = B;
  parameters.L = L;
  parameters.J_min = J_min;

  printf("Parameters for wavelet transform :\n");
  int J = s2let_j_max(&parameters);
  printf("- Wavelet parameter : %i\n", B);
  printf("- Total number of wavelets : %i\n", J);
  printf("- First wavelet scale to be used : %i\n", J_min);
  printf("- Band-limit for the map : %i\n", L);

  // Input rootname for the FITS files
  char file[100];
  char params[100];
  sprintf(params, "%d%s%d%s%d", L, "_", B, "_", J_min);
  int j, bl, offset = 0;
  printf("File root = %s\n", fileroot);

  // Read band-limits and see if multiresolution was activated
  int multires_ok = 1, monores_ok = 1;
  for (j = J; j >= J_min; j--) {
    sprintf(file, "%s%s%s%s%d%s", fileroot, "_wav_", params, "_", j, ".fits");
    printf("- Infile_wav[j=%i] = %s\n", j, file);
    bl = s2let_fits_mw_read_bandlimit(file);
    printf("  Detected bandlimit bl = %i\n", bl);
    if (bl != MIN(s2let_bandlimit(j, &parameters), L))
      multires_ok = 0;
    if (bl != L)
      monores_ok = 0;
  }
  // Read the scaling function
  sprintf(file, "%s%s%s%s", fileroot, "_scal_", params, ".fits");
  printf("- Infile_scal = %s\n", file);
  bl = s2let_fits_mw_read_bandlimit(file);
  printf("  Detected bandlimit bl = %i\n", bl);
  if (bl != MIN(s2let_bandlimit(J_min - 1, &parameters), L))
    multires_ok = 0;
  if (bl != L)
    monores_ok = 0;

  // Are the parameters and the maps all consistent ?
  if (monores_ok == 0 && multires_ok == 0) {
    printf("The parameters don't match the bandlimits of the input maps");
    printf("Neither the full or the multi-resolution algorithms are detected");
    exit(-2);
  }

  // Activate full or multi-resolution
  if (multires_ok == 1) {
    multires = 1;
    printf("Multiresolution activated\n");
  } else {
    multires = 0;
    printf("Multiresolution not activated\n");
  }

  // Allocating memory for the wavelets
  double *f_wav, *f_scal;
  if (multires) {
    s2let_transform_axisym_allocate_mw_f_wav_multires_real(
        &f_wav, &f_scal, &parameters);
  } else {
    s2let_transform_axisym_allocate_mw_f_wav_real(&f_wav, &f_scal, &parameters);
  }

  // Read the wavelets
  offset = 0;
  for (j = J_min; j <= J; j++) {
    sprintf(file, "%s%s%s%s%d%s", fileroot, "_wav_", params, "_", j, ".fits");
    if (multires)
      bl = MIN(s2let_bandlimit(j, &parameters), L);
    else
      bl = L;
    s2let_fits_mw_read_map(f_wav + offset, file, bl); // Now write the map to fits file
    offset += (2 * bl - 1) * bl;                      // Go to the next wavelet
  }
  // Read the scaling function
  sprintf(file, "%s%s%s%s", fileroot, "_scal_", params, ".fits");
  if (multires)
    bl = MIN(s2let_bandlimit(J_min - 1, &parameters), L);
  else
    bl = L;
  s2let_fits_mw_read_map(f_scal, file, bl);

  printf("Performing wavelet decomposition...");
  fflush(NULL);
  double *f = (double *)calloc(L * (2 * L - 1), sizeof(double));
  if (multires) {
    s2let_transform_axisym_wav_synthesis_mw_multires_real(
        f, f_wav, f_scal, &parameters);
  } else {
    s2let_transform_axisym_wav_synthesis_mw_real(f, f_wav, f_scal, &parameters);
  }
  printf("done\n");

  // Output the reconstruction to FITS file
  char outfile[100];
  printf("Writing the reconstructed map to a FITS file\n");
  sprintf(outfile, "%s%s%s%s", fileroot, "_recon_", params, ".fits");
  printf("- Outfile = %s\n", outfile);
  remove(outfile);                        // In case the file exists
  s2let_fits_mw_write_map(outfile, f, L); // Now write the map to fits file

  printf("--------------------------------------------------\n");

  return 0;
}
