function s2let_get_dylib
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_get_dylib
;
; PURPOSE:
;   Get the location/name of the s2let dynamic library
;
; CALLING SEQUENCE:
;   location = s2let_get_dylib()
;
;----------------------------------------------------------------------

loc = GETENV('S2LET')

r1 = file_test(loc + '/lib/libs2let.so')
r2 = file_test(loc + '/lib/libs2let.dylib')

if loc eq '' or (r1 eq 0 and r2 eq 0) then begin
   print,'You must define environment variance S2LET (base directory of s2let)'
   return, 0
endif

if r1 eq 1 then soname = loc + '/lib/libs2let.so' 
if r1 eq 0 then soname = loc + '/lib/libs2let.dylib'

return, soname

end
