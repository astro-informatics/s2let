function s2let_axisym_hpx_wav_synthesis, f_wav, verbose=verbose
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_axisym_hpx_wav_synthesis
;
; PURPOSE:
;   Reconstruct a healpix map from its wavelet transform
;
; CALLING SEQUENCE:
;   f = s2let_axisym_hpx_wav_synthesis(f_wav)
;
; INPUTS:
;   f_wav - Struc containing the wavelets maps
;
; OUTPUTS:
;   f - The reconstructed Healpix map
;
; EXAMPLE:
;
;   see s2let_hpx_demo
;
;   read_fits_map, file, f
;   f_wav = s2let_axisym_hpx_wav_analysis(f, B, L, J_min)
;   f_rec = s2let_axisym_hpx_wav_synthesis(f_wav)
;   J_max = s2let_j_max(L, b)
;   mollview, f_rec, title='Reconstructed map'
;   mollview, f_wav.scal, title='Scaling map'
;   for j=0, J_max-J_min do begin
;       mollview, f_wav.(j), title='Wavelet map '
;            +strtrim(j+1,2)+' on '+strtrim(J_max-J_min+1,2)
;   endfor
;   
; COMMENTS:
;   B, L and J_min of the input wavelets are automatically detected 
;
;----------------------------------------------------------------------


if s2let_dylib_exists() eq 0 then return, reterror('Dynamic library not found') else begin

  soname = s2let_get_dylib()

   ;nside = round(sqrt( (size(f_wav))(1) / 12.0 ))
   L = f_wav.L;
   nside = f_wav.nside
   B = f_wav.B;
   J_min = f_wav.J_min;
   J_max = f_wav.J_max;
   wavtype = f_wav.wavtype
   npix = nside2npix(nside)
   
   ;J_max = s2let_j_max(L, B)
   s2let_valid_wav_parameters, B, L, J_min, wavtype
   if keyword_set(verbose) then begin
   print, '=========================================='
   print, 's2let_axisym_hpx_wav_synthesis_real'
   print, '------------------------------------------'
   help, nside, L, B, J_min, J_max, wavtype
   endif

   f_scal = dblarr(npix)
   f_scal(0:npix-1) = f_wav.scal ;(0:npix-1, 0)
   f_wav_vec = dblarr((J_max+1-J_min)*npix)
   for j=0, J_max-J_min do begin
       f_wav_vec( j*npix : (j+1)*npix-1 ) = f_wav.(j) ;(0:npix-1, j+1)
   endfor

   f = dblarr(nside2npix(nside))
   r = call_external(soname, 's2let_idl_axisym_hpx_wav_synthesis_real', f, f_wav_vec, f_scal, nside, B, L, J_min, wavtype, /CDECL)

   if keyword_set(verbose) then print, '=========================================='

   return, f

endelse

end
