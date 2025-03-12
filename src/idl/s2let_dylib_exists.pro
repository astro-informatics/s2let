function s2let_dylib_exists
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_dylib_exists
;
; PURPOSE:
;   Check if the s2let dynamic library exists somewhere in the path
;
; CALLING SEQUENCE:
;   status = s2let_dylib_exists()
;
; OUTPUT
;   1 if found, 0 otherwise
;
;----------------------------------------------------------------------

loc = GETENV('S2LET')

r1 = file_test(loc + '/lib/libs2let.so')
r2 = file_test(loc + '/lib/libs2let.dylib')

if loc eq '' then begin
   print,'You must define environment variable S2LET (base directory of s2let)'
   return, 0
endif

if  (r1 eq 0 and r2 eq 0) then begin
   print, 'You must build the dynamic library for s2let
   return, 0
endif

return, 1
   
end
