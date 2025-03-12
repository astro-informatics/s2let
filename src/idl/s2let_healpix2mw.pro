function s2let_healpix2mw, hpxmap, lmax=lmax
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_healpix2mw
;
; PURPOSE:
;   Convert a real Healpix map to an MW map 
;   through the spherical harmonic transform 
;
; CALLING SEQUENCE:
;   mwmap = s2let_healpix2mw(hpxmap, lmax=lmax)
;
; INPUT:
;   hpxmap - The Healpix map (npix = 12*nside*nside)
;   lmax - The bandlimit of the SHA transform, thus the
;       resolution of the MW map
;
; OUTPUT
;   mwmap - The MW map (npix = lmax*(2*lmax-1))
;
;----------------------------------------------------------------------

;print, 'MAP -> ALM'
flm = s2let_hpx_map2alm_real(hpxmap, lmax)

if not keyword_set(lmax) then begin
   nside = fix(sqrt(((size(hpxmap))(1))/12.0))
   lmax = 2*nside
endif

;print, 'ALM -> MAP'
mwmap = s2let_mw_alm2map_real(flm)

return, mwmap
end
