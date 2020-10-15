// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

/*!
 * \file s2let_about.c
 * Print information about the S2LET package, including version
 * and build numbers.
 *
 * Usage: s2let_about
 *
 */

#include <stdio.h>

int main(int argc, char *argv[]) {

  printf("%s\n", "==========================================================");
  printf("%s\n", "  S2LET package");
  printf("%s\n", "  Fast Wavelets on the Sphere");
  printf("%s\n", "  By Boris Leistedt & Jason McEwen");

  printf("%s\n", "  See LICENSE.txt for license details.");

  printf("%s%s\n", "  Version: ", S2LET_VERSION);
  printf("%s%s\n", "  Build: ", S2LET_BUILD);
  printf("%s\n", "==========================================================");

  return 0;
}
