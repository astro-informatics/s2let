% s2let_demo2
% Compute and plot the harmonic tiling and the wavelet kernels.
%
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

B = 3;
zoomfactor = 1.4;
J_min = 2;
L = 128;
J = s2let_jmax(L, B);
type = 'colour';
lighting = true;

ns = ceil(sqrt(2+J-J_min+1)) ;
nx = ns - 1 + rem(2+J-J_min + 1, ns) ;
ny = ns;

[kappa kappa0] = s2let_transform_axisym_tiling(B, L, J_min);

s2let_plot_axisym_tiling(B, L, J_min);

figure('Position',[100 100 1100 700])

h=subplot(nx, ny, 1);
flm = zeros(L^2,1);
for l = 0:L-1
    flm(l^2+l+1,1) = kappa0(l+1);
end
f = ssht_inverse(flm, L, 'Reality', true);
ssht_plot_sphere(f, L, 'Type', type, 'Lighting', lighting);
zoom(zoomfactor)
v = caxis;
temp = max(abs(v));
caxis([-temp temp])
%title('k0')
locate = get(h,'title');
pos = get(locate,'position');
pos(1,2) = pos(1,2)+0.7;
pos(1,1) = pos(1,1)-0.7;
set(locate,'pos',pos);

for j = J_min:5
   h=subplot(nx, ny, j-J_min+2);
   flm = zeros(L^2,1);
    for l = 0:L-1
        flm(l^2+l+1,1) = kappa(j+1,l+1);
    end
   f = ssht_inverse(flm, L, 'Reality', true);
   ssht_plot_sphere(f, L, 'Type', type, 'Lighting', lighting);
   v = caxis;
   temp = max(abs(v));
   caxis([-temp temp])
   colormap(jet)
   %title(['j',int2str(j)-J_min+1])
   locate = get(h,'title');
   pos = get(locate,'position');
   pos(1,2) = pos(1,2)+0.7;
   pos(1,1) = pos(1,1)-0.7;
   set(locate,'pos',pos);
   zoom(zoomfactor)
end



