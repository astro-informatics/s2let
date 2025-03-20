function s2let_mw_map2alm, f
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_mw_map2alm_real
;
; PURPOSE:
;   Compute exact spherical harmonic transform of a complex MW map
;
; CALLING SEQUENCE:
;   alm = s2let_mw_map2alm(f)
;
; INPUTS
;   f - The input MW map (npix = L*(2*L-1), L is detected)
;
; OUTPUT:
;   alm - a complex array containing the spherical harmonic transform
;         alm_{el em} can be accessed at alm(ind) with ind = el*el + el + em
;
; COMMENT:
;   The resolution/bandlimit L is automatically detected
;
;----------------------------------------------------------------------


if s2let_dylib_exists() eq 1 then begin

   soname = s2let_get_dylib()

   sz = (size(f))(1)
   delta = sqrt(1 + 8*(sz))
   L = fix(( 1 + delta ) / 4)
   npix = long(L*(2*L-1))

   flm = dcomplex(dblarr(L*L))

   r = call_external(soname, 's2let_idl_mw_map2alm', flm, dcomplex(f), L, /CDECL)
   
   return, flm

endif

end
