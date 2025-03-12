pro s2let_mw_write_real_map, mapvec, file
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_mw_write_real_map
;
; PURPOSE:
;   Write a real MW map to a FITS file
;
; CALLING SEQUENCE:
;   s2let_mw_write_real_map, map, file
;
; INPUTS:
;   map - input MW map (npix=L*(2*L-1), L is detected)
;   file - filename for the FITS
;
;----------------------------------------------------------------------

L = (size(maparr))(1)

;mapvec = dblarr(L*(2*L-1))
;for el = 0, L-1 do begin
;   mapvec(el*(2*L-1):(el+1)*(2*L-1)-1) = maparr(el,*)
;endfor

sxaddpar, hdr, 'L', L

s = {data: mapvec}
mwrfits_chunks, s, file, hdr, /silent

end
