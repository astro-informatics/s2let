function s2let_mw_read_real_map, file
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_read_mw_real_map
;
; PURPOSE:
;   Read a real MW map from a FITS file
;
; CALLING SEQUENCE:
;   f = s2let_mw_read_real_map(file)
;
; INPUTS
;   file - filename for the FITS
;
; OUTPUT:
;   map - MW map (npix=L*(2*L-1), L is detected)
;
;----------------------------------------------------------------------

r = mrdfits(file, 1)

mapvec = r.(0)
sz = (size(mapvec))(1)

delta = sqrt(1 + 8*(sz))
L = ( 1 + delta ) / 4

;maparr = dblarr(L, 2*L-1)
;for el = 0, L-1 do begin
;   maparr(el,*) = mapvec(el*(2*L-1):(el+1)*(2*L-1)-1)
;endfor

return, mapvec
end
