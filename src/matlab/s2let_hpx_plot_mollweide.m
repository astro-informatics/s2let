function s2let_hpx_plot_mollweide(f)

% s2let_hpx_plot_mollweide 
% Plot a real Healpix map using Mollweide projection.
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

sz = size(f);
nside = sqrt(max(sz)/12);

[thetas, phis] = s2let_hpx_sampling_ring(nside);
          
[x, y] = ssht_mollweide(thetas, phis);

gridDelaunay = delaunay(x,y);
h = trisurf(gridDelaunay,x,y,f*0.0,f);

set(h, 'LineStyle', 'none')
axis equal
axis off
campos([0 0 1])
camup([0 1 0])


function [x, y] = ssht_mollweide(thetas, phis)
% ssht_mollweide - Compute Mollweide projection
%
% Compute Mollweide projection of spherical coordinates.
%
% Usage is given by
%
%   [x,y] = ssht_mollweide(thetas, phis)
%
% where thetas and phis are spherical coordinates and x and y are the
% projected Mollweide coordinates.
%
% Author: Jason McEwen (www.jasonmcewen.org)

MAX_ITERATIONS = 1e5;
TOL = 1e-10;

% Convert theta to longitude.
thetas = pi/2 - thetas;
phis = phis - pi;

t = thetas;
for it = 1:MAX_ITERATIONS

   dt = (t + sin(t) - pi.*sin(thetas)) ./ (1 + cos(t));
   t = t - dt;
   
   if(max(abs(dt)) < TOL)
      break;
   end
   
end
t = t/2;
x = 2 .* sqrt(2) ./ pi .* phis .* cos(t);
y = sqrt(2) .* sin(t);
