function s2let_get_mw_bandlimit, f
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_get_mw_bandlimit
;
; PURPOSE:
;   Detect the resolution/bandlimit L of an input MW map
;
; CALLING SEQUENCE:
;   L = s2let_get_mw_bandlimit(f)
;
; INPUT:
;   f - array/map of size L*(2*L-1)
;
; OUTPUT
;   L - the resolution/bandlimit
;
;----------------------------------------------------------------------


delta = sqrt(1 + 8*((size(f))(1)))
L = fix(( 1 + delta ) / 4)

return, L
end
