// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_TILING
#define S2LET_TILING
#include <ssht/ssht.h>

#define PI 3.141592653589793238462643383279502884197

#ifdef __cplusplus
extern "C" {
#endif
void s2let_switch_wavtype(int typenum);

int s2let_bandlimit(int j, const s2let_parameters_t* parameters);

int s2let_L0(int j, const s2let_parameters_t* parameters);

int s2let_j_max(const s2let_parameters_t* parameters);

void s2let_tiling_axisym_allocate(double** kappa, double** kappa0, const s2let_parameters_t* parameters);

void s2let_tiling_axisym(double* kappa, double* kappa0, const s2let_parameters_t* parameters);

void s2let_tiling_direction_allocate(S2LET_COMPLEX(double) * *s_elm, const s2let_parameters_t* parameters);

void s2let_tiling_direction(S2LET_COMPLEX(double) * s_elm, const s2let_parameters_t* parameters);

void s2let_tiling_wavelet_allocate(S2LET_COMPLEX(double) * *psi, double** phi, const s2let_parameters_t* parameters);

void s2let_tiling_wavelet(S2LET_COMPLEX(double) * psi, double* phi, const s2let_parameters_t* parameters);

double s2let_tiling_axisym_check_identity(double* kappa, double* kappa0, const s2let_parameters_t* parameters);

double s2let_tiling_direction_check_identity(S2LET_COMPLEX(double) * s_elm, const s2let_parameters_t* parameters);

double s2let_tiling_wavelet_check_identity(S2LET_COMPLEX(double) * psi, double* phi, const s2let_parameters_t* parameters);

#ifdef __cplusplus
}
#endif
#endif
