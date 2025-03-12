function s2let_hpx_map2alm_real, f, lmax
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_mw_map2alm_real
;
; PURPOSE:
;   Compute the spherical harmonic transform of a real Healpix map
;   (nb: not an exact transform due to healpix)
;
; CALLING SEQUENCE:
;   alm = s2let_hpx_map2alm_real(f, lmax)
;
; INPUTS
;   f - The input healpix map (npix = 12*nside*nside, nside is detected)
;
; OUTPUT:
;   alm - a complex array containing the spherical harmonic transform
;         alm_{el em} can be accessed at alm(ind) with ind = el*el + el + em
;
; COMMENT:
;   The resolution nside is automatically detected
;
;----------------------------------------------------------------------


if s2let_dylib_exists() eq 1 then begin

   soname = s2let_get_dylib()
   
   nside = fix(sqrt(((size(f))(1))/12.0))
print, 'nside = ', nside
print, 'lmax = ', lmax
   flm = dcomplex(dblarr(lmax^2.0))

   r = call_external(soname, 's2let_idl_hpx_map2alm_real', flm, double(f), nside, lmax, /CDECL)
   
   return, flm

endif

end
