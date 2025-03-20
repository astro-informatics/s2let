function s2let_mw_pixel_edges, L, i
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_mw_pixel_edges
;
; PURPOSE:
;   Get (theta, phi) coordinates for the corners of the i-th pixel in
;   the MW sampling
;
; CALLING SEQUENCE:
;   arr = s2let_mw_pixel_edges(L, i)
;
; INPUTS
;   L - resolution/bandlimit of the MW map
;   i - pixel index
;
; OUTPUT:
;   arr - an array containing four numbers: theta1, theta2, phi1 and phi2
;       which are the locations of the corners of the pixel
;
;----------------------------------------------------------------------

t = i / (2*L-1) 
p = i mod (2*L-1)

theta_mid = (2.0*t + 1.0) * !pi / (2.0*L - 1.0)
phi_mid = 2.0 * p * !pi / (2.0*L - 1.0)

theta1 = theta_mid - !pi / (2.0*L - 1.0)
theta2 = theta_mid + !pi / (2.0*L - 1.0)
phi1 = phi_mid - !pi / (2.0*L - 1.0)
phi2 = phi_mid + !pi / (2.0*L - 1.0)

return, [ theta1, theta2, phi1, phi2 ]

end
