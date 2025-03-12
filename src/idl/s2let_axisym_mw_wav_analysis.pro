function s2let_axisym_mw_wav_analysis, f, B, J_min, wavtype=wavtype, verbose=verbose
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_axisym_mw_wav_analysis
;
; PURPOSE:
;   Compute full resolution exact wavelet transform of a complex MW map
;
; CALLING SEQUENCE:
;   f_wav = s2let_axisym_mw_wav_analysis(f, B, J_min)
;
; INPUTS:
;   f     - input MW map (number of pixels is L*(2*L-1))
;   B     - Wavelet parameter
;   J_min - First wavelet scale to be used
;   wavtype - Wavelet type (1: scale-discretised, 2:needlets, 3: cubic splines)
;
; OUTPUTS:
;   f_wav - Struc containing the wavelets
;
; EXAMPLE:
;
;   see s2let_demo and s2let_test
;
;   f = s2let_read_mw_real_map(file)
;   if multires eq 0 then f_wav = s2let_axisym_mw_wav_analysis_real(f, B, J_min) 
;      else f_wav = s2let_axisym_mw_wav_analysis_multires_real(f, B, J_min)
;   if multires eq 0 then f_rec = s2let_axisym_mw_wav_synthesis_real(f_wav) 
;      else f_rec = s2let_axisym_mw_wav_synthesis_multires_real(f_wav)
;   L = s2let_get_mw_bandlimit(f)
;   J_max = s2let_j_max(L, B)
;   s2let_plot_mollweide, f_rec, title='Band-limited map'
;   s2let_plot_mollweide, f_wav.scal, title='Scaling map'
;   for j=0, J_max-J_min do begin
;      s2let_plot_mollweide, f_wav.(j), title='Wavelet map '
;      +strtrim(j+1,2)+' on '+strtrim(J_max-J_min+1,2)
;   endfor
;   
; COMMENTS:
;   The resolution/bandlimit L of the input map is automatically detected
;
;----------------------------------------------------------------------


if s2let_dylib_exists() eq 1 then begin

   soname = s2let_get_dylib()

   sz = (size(f))(1)
   delta = sqrt(1 + 8*(sz))
   L = fix(( 1 + delta ) / 4)
   npix = long(L*(2*L-1))
   if not keyword_set(wavtype) then wavtype = 1

   J_max = s2let_j_max(L, B)
   s2let_valid_wav_parameters, B, L, J_min, wavtype
   if keyword_set(verbose) then begin
   print, '=========================================='
   print, 's2let_axisym_mw_wav_analysis'
   print, '------------------------------------------'
   help, L, B, J_min, J_max, wavtype
   endif
   
   f_wav_vec = dcomplex(dblarr((J_max+1-J_min)*npix))
   f_scal = dcomplex(dblarr(npix))

   r = call_external(soname, 's2let_idl_transform_axisym_wav_analysis_mw', f_wav_vec, f_scal, dcomplex(f), B, L, J_min, wavtype, /CDECL)
   
   ;f_wav = dblarr(npix, J_max+2-J_min)
   ;f_wav(0:npix-1, 0) = f_scal(0:npix-1)
   f_wav = { scal: f_scal, B: B, L: L, J_min: J_min, J_max: J_max, multires: 0, maptype: 'mw', wavtype: wavtype }
   for j = J_max-J_min, 0, -1 do begin
      f_wav = CREATE_STRUCT( 'j'+strtrim(j,2), f_wav_vec( j*npix : (j+1)*npix-1 ), f_wav )
   ;   f_wav(0:npix-1, j+1) = f_wav_vec( j*npix : (j+1)*npix-1 )
   endfor

   if keyword_set(verbose) then print, '=========================================='

   return, f_wav

endif

end
