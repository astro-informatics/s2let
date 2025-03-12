function s2let_hpx_alm2map_real, flm, nside
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_hpx_alm2map_real
;
; PURPOSE:
;   Reconstruct a real Healpix map from its spherical harmonic
;   transform (nb: not an exact transform due to healpix)
;
; CALLING SEQUENCE:
;   f = s2let_hpx_alm2map_real(alm, nside)
;
; INPUT:
;   alm - a complex array containing the spherical harmonic transform
;         alm_{el em} can be accessed at alm(ind) with ind = el*el + el + em
;   nside - resolution for the output map
;
; OUTPUT
;   f - The Healpix map (npix = 12*nside*nside)
;
;----------------------------------------------------------------------


if s2let_dylib_exists() eq 1 then begin

   soname = s2let_get_dylib()
   sz = (size(flm))(1)
   L = fix(sqrt(sz))

   f = dblarr(12*nside*nside)

   r = call_external(soname, 's2let_idl_hpx_alm2map_real', f, flm, nside, L, /CDECL)
   
   return, f

endif

end
