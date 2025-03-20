function s2let_axisym_mw_wav_analysis_multires_real, f, B, J_min, wavtype=wavtype, verbose=verbose
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_axisym_mw_wav_analysis_multires_real
;
; PURPOSE:
;   Compute multiresolution exact wavelet transform of a real MW map
;
; CALLING SEQUENCE:
;   f_wav = s2let_axisym_mw_wav_analysis_multires_real(f, B, J_min)
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
   if not keyword_set(wavtype) then wavtype = 1

   J_max = s2let_j_max(L, B)
   s2let_valid_wav_parameters, B, L, J_min, wavtype
   if keyword_set(verbose) then begin
   print, '=========================================='
   print, 's2let_axisym_mw_wav_analysis_multires_real'
   print, '------------------------------------------'
   help, L, B, J_min, J_max, wavtype
   endif

   f_wav_totalsize = 0;
   for j = J_min, J_max do begin
        bl = min([s2let_get_wav_bandlimit(B, j), L])
        f_wav_totalsize = f_wav_totalsize + bl * (2 * bl - 1)
   endfor
   f_wav_vec = dblarr(f_wav_totalsize)
   bl = min([s2let_get_wav_bandlimit(B, J_min-1), L])
   f_scal = dblarr(bl*(2*bl-1))

   r = call_external(soname, 's2let_idl_axisym_mw_wav_analysis_multires_real', f_wav_vec, f_scal, double(f), B, L, J_min, wavtype, /CDECL)
   
   f_wav = { scal: f_scal, B: B, L: L, J_min: J_min, J_max: J_max, multires: 1, maptype: 'mw', wavtype: wavtype }
   offset = f_wav_totalsize
   for j = J_max, J_min, -1 do begin
      bl = min([s2let_get_wav_bandlimit(B, j), L])
      f_wav = CREATE_STRUCT( 'j'+strtrim(j,2), f_wav_vec( offset-bl*(2*bl-1) : offset-1 ), f_wav )
      offset = offset - bl*(2*bl-1)
   endfor

   if keyword_set(verbose) then print, '=========================================='

   return, f_wav

endif

end
