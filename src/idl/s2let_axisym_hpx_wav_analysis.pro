function s2let_axisym_hpx_wav_analysis, f, B, L, J_min, wavtype=wavtype, verbose=verbose
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_axisym_hpx_wav_analysis
;
; PURPOSE:
;   Compute wavelet transform of a Healpix map
;
; CALLING SEQUENCE:
;   f_wav = s2let_axisym_hpx_wav_analysis(f, B, L, J_min)
;
; INPUTS:
;   f     - input healpix nside (number of pixels is 12*nside^2)
;   B     - Wavelet parameter
;   L     - Band-limit to be used for the spherical harmonic transforms
;   J_min - First wavelet scale to be used
;   wavtype - Wavelet type (1: scale-discretised, 2:needlets, 3: cubic splines)
;
; OUTPUTS:
;   f_wav - Struc containing the wavelets
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
;   Nside of the input map is automatically detected
;
;----------------------------------------------------------------------


if s2let_dylib_exists() eq 0 then return, reterror('Dynamic library not found') else begin

   soname = s2let_get_dylib()

   nside = round(sqrt( (size(f))(1) / 12.0 ))
   npix = nside2npix(nside)
   if not keyword_set(wavtype) then wavtype = 1

   s2let_valid_wav_parameters, B, L, J_min, wavtype
   J_max = s2let_j_max(L, B)
   if keyword_set(verbose) then begin
   print, '=========================================='
   print, 's2let_axisym_hpx_wav_analysis_real'
   print, '------------------------------------------'
   help, nside, L, B, J_min, J_max, wavtype
   endif

   f_wav_vec = dblarr((J_max+1-J_min)*npix)
   f_scal = dblarr(npix)

   r = call_external(soname, 's2let_idl_axisym_hpx_wav_analysis_real', f_wav_vec, f_scal, f, nside, B, L, J_min, wavtype, /CDECL)
   
   ;f_wav = dblarr(npix, J_max+2-J_min)
   ;f_wav(0:npix-1, 0) = f_scal(0:npix-1)
   f_wav = { scal: f_scal, B: B, L: L, J_min: J_min, J_max: J_max, nside: nside, maptype: 'hpx', wavtype: wavtype }
   for j=J_max-J_min, 0, -1 do begin
      f_wav = CREATE_STRUCT( 'j'+strtrim(j,2), f_wav_vec( j*npix : (j+1)*npix-1 ), f_wav )
      ;f_wav(0:npix-1, j+1) = f_wav_vec( j*npix : (j+1)*npix-1 )
   endfor

   if keyword_set(verbose) then print, '=========================================='

   return, f_wav

endelse

end
