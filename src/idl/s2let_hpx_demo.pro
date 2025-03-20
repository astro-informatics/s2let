pro s2let_hpx_demo, B=B, L=L, J_min=J_min, wavtype=wavtype
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_hpx_demo
;
; PURPOSE:
;   Demo : test all HPX spherical harmonics and wavelet transforms
;   (not exact due to healpix) for a simulated CMB map and plot the wavelet maps
;
; OPTIONAL KEYWORDS:
;   L - The bandlimit for the spherical harmonic transforms
;   B - The wavelet parameter for the test
;   J_min - The first wavelet scale to be used for the transform
;   wavtype - Wavelet type (1: scale-discretised, 2:needlets, 3: cubic splines)
;   DEFAULT VALUES: L=192, B=7, J_min=2, wavtype=1
;
;----------------------------------------------------------------------

if not keyword_set(B) then B = 7
if not keyword_set(L) then L = 192
if not keyword_set(J_min) then J_min = 2
if not keyword_set(wavtype) then wavtype = 1

if s2let_dylib_exists() eq 1 then begin

   loc = GETENV('S2LET')
   file = loc + '/data/somecmbsimu_hpx_128.fits'
   read_fits_map, file, f

   f_wav = s2let_axisym_hpx_wav_analysis(f, B, L, J_min, wavtype=wavtype)
   f_rec = s2let_axisym_hpx_wav_synthesis(f_wav)

   J_max = s2let_j_max(L, B)
   mollview, f_rec, title='Band-limited map'
   mollview, f_wav.scal, title='Scaling map'
   for j=0, J_max-J_min do begin
      mollview, f_wav.(j), title='Wavelet map '+strtrim(j+1,2)+' on '+strtrim(J_max-J_min+1,2)
   endfor

endif

end
