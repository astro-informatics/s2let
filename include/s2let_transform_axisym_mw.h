// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_AXISYM_MW
#define S2LET_AXISYM_MW
#include <ssht/ssht.h>

#ifdef __cplusplus
extern "C" {
#endif

void s2let_transform_axisym_allocate_mw_f_wav(S2LET_COMPLEX(double) * *f_wav, S2LET_COMPLEX(double) * *f_scal, const s2let_parameters_t* parameters);
void s2let_transform_axisym_allocate_mw_f_wav_multires(S2LET_COMPLEX(double) * *f_wav, S2LET_COMPLEX(double) * *f_scal, const s2let_parameters_t* parameters);
void s2let_transform_axisym_allocate_mw_f_wav_real(double** f_wav, double** f_scal, const s2let_parameters_t* parameters);
void s2let_transform_axisym_allocate_mw_f_wav_multires_real(double** f_wav, double** f_scal, const s2let_parameters_t* parameters);

void s2let_transform_axisym_wav_analysis_mw(S2LET_COMPLEX(double) * f_wav, S2LET_COMPLEX(double) * f_scal, const S2LET_COMPLEX(double) * f, const s2let_parameters_t* parameters);
void s2let_transform_axisym_wav_analysis_adjoint_mw(S2LET_COMPLEX(double) * f, const S2LET_COMPLEX(double) * f_wav, const S2LET_COMPLEX(double) * f_scal, const s2let_parameters_t* parameters);
void s2let_transform_axisym_wav_synthesis_mw(S2LET_COMPLEX(double) * f, const S2LET_COMPLEX(double) * f_wav, const S2LET_COMPLEX(double) * f_scal, const s2let_parameters_t* parameters);
void s2let_transform_axisym_wav_synthesis_adjoint_mw(S2LET_COMPLEX(double) * f_wav, S2LET_COMPLEX(double) * f_scal, const S2LET_COMPLEX(double) * f, const s2let_parameters_t* parameters);

void s2let_transform_axisym_wav_analysis_mw_multires(S2LET_COMPLEX(double) * f_wav, S2LET_COMPLEX(double) * f_scal, const S2LET_COMPLEX(double) * f, const s2let_parameters_t* parameters);
void s2let_transform_axisym_wav_analysis_adjoint_mw_multires(S2LET_COMPLEX(double) * f, const S2LET_COMPLEX(double) * f_wav, const S2LET_COMPLEX(double) * f_scal, const s2let_parameters_t* parameters);
void s2let_transform_axisym_wav_synthesis_mw_multires(S2LET_COMPLEX(double) * f, const S2LET_COMPLEX(double) * f_wav, const S2LET_COMPLEX(double) * f_scal, const s2let_parameters_t* parameters);
void s2let_transform_axisym_wav_synthesis_adjoint_mw_multires(S2LET_COMPLEX(double) * f_wav, S2LET_COMPLEX(double) * f_scal, const S2LET_COMPLEX(double) * f, const s2let_parameters_t* parameters);

void s2let_transform_axisym_wav_analysis_mw_real(double* f_wav, double* f_scal, const double* f, const s2let_parameters_t* parameters);
void s2let_transform_axisym_wav_analysis_adjoint_mw_real(double* f, const double* f_wav, const double* f_scal, const s2let_parameters_t* parameters);
void s2let_transform_axisym_wav_synthesis_mw_real(double* f, const double* f_wav, const double* f_scal, const s2let_parameters_t* parameters);
void s2let_transform_axisym_wav_synthesis_adjoint_mw_real(double* f_wav, double* f_scal, const double* f, const s2let_parameters_t* parameters);

void s2let_transform_axisym_wav_analysis_mw_multires_real(double* f_wav, double* f_scal, const double* f, const s2let_parameters_t* parameters);
void s2let_transform_axisym_wav_analysis_adjoint_mw_multires_real(double* f, const double* f_wav, const double* f_scal, const s2let_parameters_t* parameters);
void s2let_transform_axisym_wav_synthesis_mw_multires_real(double* f, const double* f_wav, const double* f_scal, const s2let_parameters_t* parameters);
void s2let_transform_axisym_wav_synthesis_adjoint_mw_multires_real(double* f_wav, double* f_scal, const double* f, const s2let_parameters_t* parameters);

void s2let_transform_axisym_wav_hardthreshold_real(double* g_wav, const double* threshold, const s2let_parameters_t* parameters);
void s2let_transform_axisym_wav_hardthreshold_multires_real(double* g_wav, const double* threshold, const s2let_parameters_t* parameters);

#ifdef __cplusplus
}
#endif
#endif
