pro s2let_test, L=L, B=B, J_min=J_min, multires=multires, wavtype=wavtype
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_test
;
; PURPOSE:
;   Test the exactness of all MW spherical harmonics and wavelet
;   transform for real and complex signals on the sphere
;
; OPTIONAL KEYWORDS:
;   L - The bandlimit/resolution for the test
;   B - The wavelet parameter for the test
;   J_min - The first wavelet scale to be used for the transform
;   multires - multiresolution flag (1 if yes, else 0)
;   wavtype - Wavelet type (1: scale-discretised, 2:needlets, 3: cubic splines)
;   DEFAULT VALUES: L=64, B=3, J_min=2, multires=0, wavtype=1
;
;----------------------------------------------------------------------

if not keyword_set(L) then L = 64
if not keyword_set(multires) then multires = 0
if not keyword_set(B) then B = 3
if not keyword_set(J_min) then J_min = 2
if not keyword_set(wavtype) then wavtype = 1

if s2let_dylib_exists() eq 1 then begin

   loc = GETENV('S2LET')

   flm = (2.0 * dcomplex(randomn(1,L*L), randomn(1,L*L))- 1.0)
   for el = 0, L-1 do begin
      ind = el*el + el
      flm(ind) = REAL_PART(flm(ind))
      for em = 1, el do begin
         ind = el*el + el + em
         indb = el*el + el - em
         flm(indb) = (-1.0)^em * conj( flm(ind) )
      endfor
   endfor

   f = s2let_mw_alm2map_real(flm)

   if multires eq 0 then f_wav = s2let_axisym_mw_wav_analysis_real(f, B, J_min, wavtype=wavtype) else f_wav = s2let_axisym_mw_wav_analysis_multires_real(f, B, J_min, wavtype=wavtype)
   if multires eq 0 then f_rec = s2let_axisym_mw_wav_synthesis_real(f_wav) else f_rec = s2let_axisym_mw_wav_synthesis_multires_real(f_wav)

   flm_rec = s2let_mw_map2alm_real(f_rec)

   print, 'Maximum absolute error on real transform :', max(abs(flm_rec-flm))

   flm2 = (2.0 * dcomplex(randomn(1,L*L), randomn(1,L*L))- 1.0)

   f2 = s2let_mw_alm2map(flm2)

   if multires eq 0 then f_wav2 = s2let_axisym_mw_wav_analysis(f2, B, J_min, wavtype=wavtype) else f_wav = s2let_axisym_mw_wav_analysis_multires(f2, B, J_min, wavtype=wavtype)

   if multires eq 0 then f_rec2 = s2let_axisym_mw_wav_synthesis(f_wav2) else f_rec2 = s2let_axisym_mw_wav_synthesis_multires(f_wav2)

   flm_rec2 = s2let_mw_map2alm(f_rec2)

   print, 'Maximum absolute error on complex transform :', max(abs(flm_rec2-flm2))

endif

end
