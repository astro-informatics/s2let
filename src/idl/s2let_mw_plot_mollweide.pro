pro s2let_mw_plot_mollweide, maporfile, nlevels=nlevels, title=title, charsize=charsize
;+
; S2LET package - Copyright (C) 2012 
; Boris Leistedt & Jason McEwen
;
; NAME:
;   s2let_mw_plot_mollweide
;
; PURPOSE:
;   Plot a real MW map (from data of FITS file) under Mollweide projection
;
; CALLING SEQUENCE:
;   s2let_mw_plot_mollweide, mapfile, title='My map'
;
; INPUTS
;   maporfile - filename for the FITS or map read with s2let_read_mw_real_map
;
;----------------------------------------------------------------------


if not keyword_set(title) then title='MW map'
if not keyword_set(charsize) then charsize=2.0
if not keyword_set(nlevels) then nlevels=100
if total(valid_num(maporfile)) eq 0 then map = s2let_mw_read_real_map(maporfile) else map = maporfile


thecolors = 18 + indgen(255-18)
ncolors = n_elements(thecolors)
h = (max(map) - min(map))/double(ncolors-1)
valrange = min(map) + indgen(ncolors)*h

sz = (size(map))(1)
delta = sqrt(1 + 8*(sz))
L = long(fix(( 1 + delta ) / 4))

maparr = dblarr(L, 2*L-1)
for el = 0, L-1 do begin
   maparr(el, *) = map(el*(2*L-1):(el+1)*(2*L-1)-1)
endfor

s2let_mw_sampling, L, thetas, phis

MAX_ITERATIONS = 1e2;
TOL = 1e-5;
; Convert theta to longitude.
;thetas = !pi/2 - thetas;;
;phis = phis - !pi;
;t = thetas;
;for it = 1, MAX_ITERATIONS do begin
;   dt = (t + sin(t) - !pi*sin(thetas)) / (1 + cos(t));
;   t = t - dt;
;   if max(abs(dt)) lt TOL then break
;endfor
;t = t/2.0;
;x = 2 * sqrt(2) / !pi * cos(t) # phis;
;y = sqrt(2) * sin(t) # (phis/phis);

;maparr_b = rebin(maparr, L*10, (2*L-1)*10, /sample)
;x_b = rebin(x, L*10, (2*L-1)*10, /sample)
;y_b = rebin(y, L*10, (2*L-1)*10, /sample)

;tvimage, reverse(transpose(maparr_b),2), transpose(x_b),
;transpose(y_b), /KEEP_ASPECT_RATIO, /nointerpolation

plot, [-2.9,2.9], [-1.42,1.42], xstyle=5, ystyle=5, title=title, color=17, charsize=charsize, /nodata, xrange=[-2.9,2.9], yrange=[-1.42,1.42]

cols = floor(10*(RANDOMN(1, L*(2*L-1))))
npixtot =long(L)*(2*long(L)-1)
i = long(0)
repeat begin

   arr = s2let_mw_pixel_edges(L, i) ; [t1, t2, p1, p2]
   thetas = [arr(0), arr(1)]  
   phis = [arr(2), arr(3)]
   thetas = !pi/2.0 - thetas      
   phis = phis - !pi  
   t = thetas
   for it = 1, MAX_ITERATIONS do begin
      dt = (t + sin(t) - !pi*sin(thetas)) / (1 + cos(t)) 
      t = t - dt                                         
      if max(abs(dt)) lt TOL then break
   endfor
   t = t/2.0                          
   x = 2 * sqrt(2) / !pi * cos(t) # phis 
   y = sqrt(2) * sin(t) # (0.0*phis+1.0);

   x = [x(0:1), x(3), x(2)]
   y = [y(0:1), y(3), y(2)]

   coli = (where( abs(valrange-map(i)) eq (min(abs(valrange-map(i))))(0) ))(0)
   
   polyfill, x, y, color=thecolors(coli)
   
   i = i + 1
endrep until i eq npixtot

;contour, maparr, x, y,  /fill, nlevels=nlevels, xstyle=4, ystyle=4, title=title, color=17, charsize=charsize, /nodata

end
