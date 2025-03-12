pro s2let_valid_wav_parameters, B, L, J_min, wavtype
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_valid_wav_parameters
;
; PURPOSE:
;   Check and test a set of wavelet parameters
;
; CALLING SEQUENCE:
;   status = s2let_valid_wav_parameters(B, L, J_min)
;
; INPUTS:
;   B     - Wavelet parameter
;   L     - Band-limit to be used for the spherical harmonic transforms
;   J_min - First wavelet scale to be used
;   wavtype - Wavelet type (1: scale-discretised, 2:needlets, 3: cubic splines)
;
;----------------------------------------------------------------------

if not valid_num(B,/integer) or B lt 2 then stop, 'Error: Parameter B must be a positive integer'
if not valid_num(L,/integer) or L lt 2 then stop, 'Error: Parameter L must be a positive integer'
if not valid_num(J_min,/integer) or J_min lt 0 then stop, 'Error: Parameter J_min must be a positive integer'
if not valid_num(wavtype,/integer) or wavtype gt 3 or wavtype lt 1 then stop, 'Error: Parameter wavtype must be 1, 2 or 3'
J_max = s2let_j_max(L, B)
msg = 'Error: Parameter B must be greater than J_max='+strtrim(J_max,2)
if J_max lt 1 then stop, msg
msg = 'Error: Parameter J_min must be lower than J_max='+strtrim(J_max,2)
if J_min ge J_max then stop, msg

end
