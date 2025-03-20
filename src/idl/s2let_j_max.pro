function s2let_j_max, L, B
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_j_max
;
; PURPOSE:
;   Compute the maximum wavelet to be used, given a bandlimit L and a
;   wavelet parameter B
;
; CALLING SEQUENCE:
;   J_max = s2let_j_max(L, B)
;
; OUTPUT
;   Compute the maximum wavelet to be used, given a bandlimit L and a
;   wavelet parameter B;
;
;----------------------------------------------------------------------


return, fix(ceil(alog10(L) / alog10(B)))

end
