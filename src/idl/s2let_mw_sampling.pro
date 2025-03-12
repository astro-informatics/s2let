pro s2let_mw_sampling, L, thetas, phis
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_mw_sampling
;
; PURPOSE:
;   Compute the coordinates of the nodes of the MW sampling for a
;   given bandlimit/resolution
;
; CALLING SEQUENCE:
;   s2let_mw_sampling, L, thetas, phis
;
; INPUT:
;   L - resolution/bandlimit
;
; OUTPUT
;   thetas - The theta coordinates of the nodes
;   phis - The phi coordinates of the nodes
;
;----------------------------------------------------------------------

thetas = dblarr(L)
phis = dblarr(2*L-1)

;t = long(0)
;repeat begin
for t = 0, L-1 do begin
   thetas(t) = (2.0*t + 1.0) * !pi / (2.0*L - 1.0)
endfor
;   t = t + 1
;endrep until t eq L - 1

;p = long(0)
;repeat begin
for p = 0, 2*L-2 do begin
   phis(p) = 2.0 * p * !pi / (2.0*L - 1.0);
endfor
;   p = p + 1
;endrep until p eq 2 * L - 2

end
