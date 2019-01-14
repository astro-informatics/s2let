// S2LET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef S2LET_IDL_MW
#define S2LET_IDL_MW

/*!
 * IDL interface to s2let_mw_axisym_wav_analysis_mw_real
 */
int s2let_idl_transform_axisym_wav_analysis_mw_real(int argc, void* argv[]);

/*!
 * IDL interface to s2let_mw_axisym_wav_synthesis_mw_real
 */
int s2let_idl_transform_axisym_wav_synthesis_mw_real(int argc, void* argv[]);

/*!
 * IDL interface to s2let_mw_axisym_wav_analysis_mw
 */
int s2let_idl_transform_axisym_wav_analysis_mw(int argc, void* argv[]);

/*!
 * IDL interface to s2let_mw_axisym_wav_synthesis_mw
 */
int s2let_idl_transform_axisym_wav_synthesis_mw(int argc, void* argv[]);


/*!
 * IDL interface to s2let_mw_axisym_wav_analysis_mw_multires_real
 */
int s2let_idl_transform_axisym_wav_analysis_mw_multires_real(int argc, void* argv[]);

/*!
 * IDL interface to s2let_mw_axisym_wav_synthesis_mw_multires_real
 */
int s2let_idl_transform_axisym_wav_synthesis_mw_multires_real(int argc, void* argv[]);

/*!
 * IDL interface to s2let_mw_axisym_wav_analysis_mw_multires
 */
int s2let_idl_transform_axisym_wav_analysis_mw_multires(int argc, void* argv[]);

/*!
 * IDL interface to s2let_mw_axisym_wav_synthesis_mw_multires
 */
int s2let_idl_transform_axisym_wav_synthesis_mw_multires(int argc, void* argv[]);

/*!
 * IDL interface to s2let_mw_map2alm_real
 */
int s2let_idl_transform_map2alm_real(int argc, void* argv[]);   

/*!
 * IDL interface to s2let_mw_alm2map_real
 */
int s2let_idl_transform_alm2map_real(int argc, void* argv[]);

/*!
 * IDL interface to s2let_mw_map2alm
 */
int s2let_idl_transform_map2alm(int argc, void* argv[]);   

/*!
 * IDL interface to s2let_mw_alm2map
 */
int s2let_idl_transform_alm2map(int argc, void* argv[]);

#endif
