function s2let_get_wav_bandlimit, B, j
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_get_wav_bandlimit
;
; PURPOSE:
;   Compute the bandlimit of the j-th wavelet constructed with parameter B
;
; CALLING SEQUENCE:
;   bl = s2let_get_wav_bandlimit(B, j)
;
; OUTPUT
;   The bandlimit of the j-th wavelet constructed with parameter B
;
;----------------------------------------------------------------------

return, ceil(B^float(j+1))
end
