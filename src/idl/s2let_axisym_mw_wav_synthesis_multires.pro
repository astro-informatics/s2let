function s2let_axisym_mw_wav_synthesis_multires, f_wav, verbose=verbose
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_axisym_mw_wav_synthesis_multires
;
; PURPOSE:
;   Reconstruct a complex MW map from its multiresolution exact wavelet transform
;
; CALLING SEQUENCE:
;   f = s2let_axisym_mw_wav_synthesis_multires(f_wav)
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
      
   s2let_valid_wav_parameters, B, L, J_min, wavtype
   if keyword_set(verbose) then begin
   print, '=========================================='
   print, 's2let_axisym_mw_wav_synthesis_multires_real'
   print, '------------------------------------------'
   help, L, B, J_min, J_max
   endif
   
   bl = min([s2let_get_wav_bandlimit(B, J_min-1), L])
   f_scal = f_wav.scal 
   f_wav_totalsize = 0           ;
   for j = J_min, J_max do begin
        bl = min([s2let_get_wav_bandlimit(B, j), L])
        f_wav_totalsize = f_wav_totalsize + bl * (2 * bl - 1)
   endfor
   f_wav_vec = dcomplex(dblarr(f_wav_totalsize))
   offset = 0
   for j = J_min, J_max do begin
       bl = min([s2let_get_wav_bandlimit(B, j), L])  
       f_wav_vec( offset : offset + bl*(2*bl-1)-1 ) = f_wav.(j-J_min)
       offset = offset + bl * (2 * bl - 1)
    endfor

   f = dcomplex(dblarr(L*(2*L-1)))
   r = call_external(soname, 's2let_idl_transform_axisym_wav_synthesis_mw_multires_real', f, dcomplex(f_wav_vec), dcomplex(f_scal), B, L, J_min, wavtype, /CDECL)

   if keyword_set(verbose) then print, '=========================================='

   return, f

endif

end
