// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#ifndef S2LET_AXISYM_MW
#define S2LET_AXISYM_MW

#include <complex.h>

void s2let_transform_axisym_allocate_mw_f_wav(complex double **f_wav, complex double **f_scal, const s2let_parameters_t *parameters);
void s2let_transform_axisym_allocate_mw_f_wav_multires(complex double **f_wav, complex double **f_scal, const s2let_parameters_t *parameters);
void s2let_transform_axisym_allocate_mw_f_wav_real(double **f_wav, double **f_scal, const s2let_parameters_t *parameters);
void s2let_transform_axisym_allocate_mw_f_wav_multires_real(double **f_wav, double **f_scal, const s2let_parameters_t *parameters);

void s2let_transform_axisym_wav_analysis_mw(complex double *f_wav, complex double *f_scal, const complex double *f, const s2let_parameters_t *parameters);
void s2let_transform_axisym_wav_analysis_adjoint_mw(complex double *f, const complex double *f_wav, const complex double *f_scal, const s2let_parameters_t *parameters);
void s2let_transform_axisym_wav_synthesis_mw(complex double *f, const complex double *f_wav, const complex double *f_scal, const s2let_parameters_t *parameters);
void s2let_transform_axisym_wav_synthesis_adjoint_mw(complex double *f_wav, complex double *f_scal, const complex double *f, const s2let_parameters_t *parameters);

void s2let_transform_axisym_wav_analysis_mw_multires(complex double *f_wav, complex double *f_scal, const complex double *f, const s2let_parameters_t *parameters);
void s2let_transform_axisym_wav_analysis_adjoint_mw_multires(complex double *f, const complex double *f_wav, const complex double *f_scal, const s2let_parameters_t *parameters);
void s2let_transform_axisym_wav_synthesis_mw_multires(complex double *f, const complex double *f_wav, const complex double *f_scal, const s2let_parameters_t *parameters);
void s2let_transform_axisym_wav_synthesis_adjoint_mw_multires(complex double *f_wav, complex double *f_scal, const complex double *f, const s2let_parameters_t *parameters);

void s2let_transform_axisym_wav_analysis_mw_real(double *f_wav, double *f_scal, const double *f, const s2let_parameters_t *parameters);
void s2let_transform_axisym_wav_analysis_adjoint_mw_real(double *f, const double *f_wav, const double *f_scal, const s2let_parameters_t *parameters);
void s2let_transform_axisym_wav_synthesis_mw_real(double *f, const double *f_wav, const double *f_scal, const s2let_parameters_t *parameters);
void s2let_transform_axisym_wav_synthesis_adjoint_mw_real(double *f_wav, double *f_scal, const double *f, const s2let_parameters_t *parameters);

void s2let_transform_axisym_wav_analysis_mw_multires_real(double *f_wav, double *f_scal, const double *f, const s2let_parameters_t *parameters);
void s2let_transform_axisym_wav_analysis_adjoint_mw_multires_real(double *f, const double *f_wav, const double *f_scal, const s2let_parameters_t *parameters);
void s2let_transform_axisym_wav_synthesis_mw_multires_real(double *f, const double *f_wav, const double *f_scal, const s2let_parameters_t *parameters);
void s2let_transform_axisym_wav_synthesis_adjoint_mw_multires_real(double *f_wav, double *f_scal, const double *f, const s2let_parameters_t *parameters);

void s2let_transform_axisym_wav_hardthreshold_real(double *g_wav, const double *threshold, const s2let_parameters_t *parameters);
void s2let_transform_axisym_wav_hardthreshold_multires_real(double *g_wav, const double *threshold, const s2let_parameters_t *parameters);

#endif
