function s2let_axisym_mw_wav_synthesis_real, f_wav;, B, J_min, verbose=verbose
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_axisym_mw_wav_synthesis_real
;
; PURPOSE:
;   Reconstruct a real MW map from its full resolution exact wavelet transform
;
; CALLING SEQUENCE:
;   f = s2let_axisym_mw_wav_synthesis_real(f_wav)
;
; INPUTS:
;   f_wav - input exact wavelet transfrom 
;
; OUTPUTS:
;   f     - Reconstructed MW map
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
;   The wavelet parameters B/J_min are automatically detected
;
;----------------------------------------------------------------------


if s2let_dylib_exists() eq 1 then begin

   soname = s2let_get_dylib()

   ;sz = (size(f_wav))(1)
   ;delta = sqrt(1 + 8*(sz))
   ;L = fix(( 1 + delta ) / 4)
   L = f_wav.L;
   B = f_wav.B;
   J_min = f_wav.J_min;
   J_max = f_wav.J_max;
   wavtype = f_wav.wavtype
   
   npix = long(L)*(2*long(L)-1)
   
   s2let_valid_wav_parameters, B, L, J_min, wavtype
   ;J_max = s2let_j_max(L, B)
   if keyword_set(verbose) then begin
   print, '=========================================='
   print, 's2let_axisym_mw_wav_synthesis_real'
   print, '------------------------------------------'
   help, L, B, J_min, J_max
   endif

   f_scal = dblarr(npix)
   f_scal(0:npix-1) = f_wav.scal ;f_wav(0:npix-1, 0)
   f_wav_vec = dblarr((J_max+1-J_min)*npix)
   for j=0, J_max-J_min do begin
       ;jnm = 'j'+strtrim(j,2)
       f_wav_vec( j*npix : (j+1)*npix-1 ) = f_wav.(j);jnm ;f_wav(0:npix-1, j+1)
   endfor

   f = dblarr(npix)
   r = call_external(soname, 's2let_idl_transform_axisym_wav_synthesis_mw_real', f, double(f_wav_vec), double(f_scal), B, L, J_min, wavtype, /CDECL)

   if keyword_set(verbose) then print, '=========================================='

   return, f

endif

end
