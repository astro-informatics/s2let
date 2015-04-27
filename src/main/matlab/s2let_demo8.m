% s2let_demo8
% -----------------------------------------------------------
% Curvelets: 
% Compute and plot the harmonic tiling and the wavelet kernels.
%
% ------ About this demo -------- 
% It has a similar structure as s2let_demo5.m
% but here s2let_curvelet_tiling(B, L, N, Spin, J_min) is called 
% instead of s2let_wavelet_tiling(B, L, N, Spin, J_min) for computing directional wavelets.
%
% Same functionality as s2let_demo7, but demo8 interfaces with C. 
% See the details of the function "s2let_curvelet_tiling" in s2let_tiling.c
% -----------------------------------------------------------
% Log: 
% -  constructed by Jennifer Y H Chan on 10th Dec 2014  
% -----------------------------------------------------------
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details
% 
% Modified S2LET package to perform Curvelets on the Sphere.
% -----------------------------------------------------------

Spin = 0;
B = 2;   % for dyadic sampling 
L = 64;
N= L;  % Since m=l, the azimuthal band limit N = overall band limit L

% minimum and maximum scale probed by wavelets 
J_min = 2;
J = s2let_jmax(L, B);  % J = ceil(log L/ log B)

zoomfactor = 1.4;
plot_caxis_scale = 2
type = 'colour';
lighting = true;

%ns = ceil(sqrt(2+J-J_min+1)) ;
%nx = ns - 1 + rem(2+J-J_min + 1, ns) ;
%ny = ns;
nx = 3
ny = 4

maxfigs = nx*ny;
pltroot = '../../../figs'
configstr = ['Spin',int2str(Spin),'_N',int2str(N),'_L',int2str(L),'_B',int2str(B),'_Jmin',int2str(J_min)]

% Compute curvelets coefficients (s2let_tiling.c):
[psi_lm phi_l] = s2let_curvelet_tiling(B, L, N, Spin, J_min);

% Plot of tiling of curvelets in harmonic space:]
% s2let_plot_cur_tiling(B, L, J_min);


% Plot curvelets on the sphere
% Rotate the sphere such that x-axis -> z-axis 
% (i.e. curvelets lie on the sphere top)
gamma = 0
beta = pi/2
alpha = pi
% Precompute Wigner small-d functions for rotation 
d = zeros(L, 2*L-1, 2*L-1);
d(1,:,:) = ssht_dl(squeeze(d(1,:,:)), L, 0, beta);
for el = 1:L-1
    d(el+1,:,:) = ssht_dl(squeeze(d(el,:,:)), L, el, beta);
end

figure('Position',[100 100 1200 600])
ind = 0;
for j = J_min:J
   flm = psi_lm(:,j+1);
   
% Rotate spherical harmonic coefficients
   flm_rot = ssht_rotate_flm(flm, d, alpha, gamma);

   if Spin == 0
% Compute the non-rotated function
%       f = ssht_inverse(flm, L, 'Reality', true);
% Compute the rotated function
       f_rot = ssht_inverse(flm_rot, L, 'Reality', true);
       ind = ind + 1;
       if ind <= maxfigs
           h = subplot(ny, nx, ind);
% Plot the non-rotated function
%           ssht_plot_sphere(f, L, 'Type', type, 'Lighting', lighting);   
% Plot the rotated function
           ssht_plot_sphere(f_rot, L, 'Type', type, 'Lighting', lighting);
           title(h, ['Wavelet j = ',int2str(j-J_min+1)])
           locate = get(h,'title');
           pos = get(locate,'position');
           pos(1,2) = pos(1,2)+0.7;
           pos(1,1) = pos(1,1)-0.7;
           set(locate,'pos',pos);
           v = caxis;
           temp = max(abs(v));
           caxis([-temp temp]*plot_caxis_scale)
           zoom(zoomfactor)

       end
   end
   if Spin > 0
% Compute the non-rotated function
%       f = ssht_inverse(flm, L, 'spin', Spin);
% Compute the rotated function
       f_rot = ssht_inverse(flm_rot, L, 'spin', Spin);
       ind = ind + 1;
       if ind <= maxfigs
           h = subplot(ny, nx, ind);
% Plot the real components of the non-rotated function:        
%           ssht_plot_sphere(real(f), L, 'Type', type, 'Lighting', lighting);    
% Plot the real components of the rotated function:  
           ssht_plot_sphere(real(f_rot), L, 'Type', type, 'Lighting', lighting);
           title(h, ['Wavelet j = ',int2str(j-J_min+1), ', real part'])
           locate = get(h,'title');
           pos = get(locate,'position');
           pos(1,2) = pos(1,2)+0.7;
           pos(1,1) = pos(1,1)-0.7;
           set(locate,'pos',pos);
           v = caxis;
           temp = max(abs(v));
           caxis([-temp temp]*plot_caxis_scale)
           zoom(zoomfactor)

       end
       ind = ind + 1;
       if ind <= maxfigs
           h = subplot(ny, nx, ind);
% Plot the imaginary components of the non-rotated function:     
%           ssht_plot_sphere(imag(f), L, 'Type', type, 'Lighting', lighting);
% Plot the imaginary components of the rotated function:   
           ssht_plot_sphere(imag(f_rot), L, 'Type', type, 'Lighting', lighting);
           title(h, ['Wavelet j = ',int2str(j-J_min+1), ', imag part'])
           locate = get(h,'title');
           pos = get(locate,'position');
           pos(1,2) = pos(1,2)+0.7;
           pos(1,1) = pos(1,1)-0.7;
           set(locate,'pos',pos);
           v = caxis;
           temp = max(abs(v));
           caxis([-temp temp]*plot_caxis_scale)
           zoom(zoomfactor)
           
       end
       
       ind = ind + 1;
       if ind <= maxfigs
           h = subplot(ny, nx, ind);
% Plot the positive non-rotated function:     
%           ssht_plot_sphere(abs(f), L, 'Type', type, 'Lighting', lighting);
% Plot the positive rotated function:     
           ssht_plot_sphere(abs(f_rot), L, 'Type', type, 'Lighting', lighting);          
           title(h, ['Wavelet j = ',int2str(j-J_min+1), ', abs part'])
           locate = get(h,'title');
           pos = get(locate,'position'); 
           pos(1,2) = pos(1,2)+0.7;
           pos(1,1) = pos(1,1)-0.7;
           set(locate,'pos',pos);
           v = caxis;
           temp = max(abs(v));
           caxis([-temp temp]*plot_caxis_scale)
           zoom(zoomfactor)
       end
       
   end
end


colormap(jet)
fname = ['s2let_demo8_', configstr, '_cur_jet.png']
print('-r200', '-dpng', fname)
%colormap(hot)
%fname = ['s2let_demo8_', configstr,'_cur_hot.png']
%print('-r200', '-dpng', fname)



% Scaling function - axisymmetric (m=0)
% Note the index convention of flm: index= el .* el + el + m + 1
figure('Position',[100 100 300 300])
h=subplot(1, 1, 1);
flm = zeros(L^2,1);
for l = 0:L-1
   % Curvelet scaling coefficients
   flm(l^2+l+1,1) = phi_l(l+1);
end
% Compute the scaling coefficients to obatin the scaling function
f = ssht_inverse(flm, L, 'Reality', true);
% Plot the scaling function on the sphere 
ssht_plot_sphere(f, L, 'Type', type, 'Lighting', lighting);
title(h,'Scaling fct')
locate = get(h,'title');
pos = get(locate,'position');
pos(1,2) = pos(1,2)+0.7;
pos(1,1) = pos(1,1)-0.7;
set(locate,'pos',pos);
zoom(1.2)
v = caxis;
temp = max(abs(v));
caxis([-temp temp]*plot_caxis_scale)


colormap(jet)
% fname = [pltroot,'/s2let_demo8_', configstr, '_scal_jet.png']
fname = ['s2let_demo8_', configstr, '_scal_jet.png']
print('-r200', '-dpng', fname)
%colormap(hot)
%% fname = [pltroot,'/s2let_demo8_', configstr, '_scal_hot.png']
%fname = ['s2let_demo8_', configstr, '_scal_hot.png']
%print('-r200', '-dpng', fname)
