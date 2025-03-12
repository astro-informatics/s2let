function s2let_mw_alm2map_real, flm
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_mw_alm2map_real
;
; PURPOSE:
;   Reconstruct a real MW map from its exact spherical harmonic transform
;
; CALLING SEQUENCE:
;   f = s2let_mw_alm2map_real(alm)
;
; INPUT:
;   alm - a complex array containing the spherical harmonic transform
;         alm_{el em} can be accessed at alm(ind) with ind = el*el + el + em
;
; OUTPUT
;   f - The MW map (npix = L*(2*L-1), L is detected)
;
; COMMENT:
;   The resolution/bandlimit L is automatically detected
;
;----------------------------------------------------------------------


if s2let_dylib_exists() eq 1 then begin

   soname = s2let_get_dylib()
   sz = (size(flm))(1)
   delta = sqrt(1 + 8*(sz))
   L = long(sqrt(sz))

   f = dblarr(L*(2*L-1))

   r = call_external(soname, 's2let_idl_mw_alm2map_real', f, flm, L, /CDECL)
   
   return, f

endif

end
