pro s2let_demo2, L=L, B=B, J_min=J_min, multires=multires, wavtype=wavtype
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_demo2
;
; PURPOSE:
;   Demo : read Healpix map, convert it to MW map, then test the
;   exactness of all MW spherical harmonics and wavelet
;   transforms and plot the wavelet maps
;
; OPTIONAL KEYWORDS:
;   L - The bandlimit for the spherical harmonic transforms
;   B - The wavelet parameter for the test
;   J_min - The first wavelet scale to be used for the transform
;   multires - multiresolution flag (1 if yes, else 0)
;   wavtype - Wavelet type (1: scale-discretised, 2:needlets, 3: cubic splines)
;   DEFAULT VALUES: L=128, B=3, J_min=2, multires=0, wavtype=1
;
;----------------------------------------------------------------------


if not keyword_set(L) then L = 394
if not keyword_set(multires) then multires = 0
if not keyword_set(B) then B = 2
if not keyword_set(J_min) then J_min = 3
if not keyword_set(wavtype) then wavtype = 1

if s2let_dylib_exists() eq 1 then begin

   loc = GETENV('S2LET')
   file = loc + '/data/somecmbsimu_hpx_128.fits'
   ;file = '/Users/bl/Dropbox/Melissa/spectre/out/lsspaperbins/dr6quasars_map_z065145_256_Q.fits'
   read_fits_map, file, f_hpx

   f = s2let_healpix2mw(f_hpx, lmax=L)

   if multires eq 0 then f_wav = s2let_axisym_mw_wav_analysis_real(f, B, J_min, wavtype=wavtype, verbose=verbose) else f_wav = s2let_axisym_mw_wav_analysis_multires_real(f, B, J_min, wavtype=wavtype, verbose=verbose)

   if multires eq 0 then f_rec = s2let_axisym_mw_wav_synthesis_real(f_wav) else f_rec = s2let_axisym_mw_wav_synthesis_multires_real(f_wav)

   sz = (size(f))(1)
   delta = sqrt(1 + 8*(sz))
   L = fix(( 1 + delta ) / 4)
   J_max = s2let_j_max(L, b)

   ns = ceil(sqrt(3+J_max-J_min)) 
   nx = ns
   if 3+J_max-J_min le ns*(ns-1) then begin
      ny = ns - 1 
      ;window, xsize=1200, ysize=533
   endif else begin 
      ny = ns
      ;window, xsize=1200, ysize=800
   endelse

   !P.MULTI=[0,nx,ny]

   s2let_mw_plot_mollweide, f_rec, title='Band-limited map'
   s2let_mw_plot_mollweide, f_wav.scal, title='Scaling map'
   for j=0, J_max-J_min do begin
      s2let_mw_plot_mollweide, f_wav.(j), title='Wavelet map '+strtrim(j+1,2)+' on '+strtrim(J_max-J_min+1,2)
   endfor

   !P.MULTI=0

endif

end
