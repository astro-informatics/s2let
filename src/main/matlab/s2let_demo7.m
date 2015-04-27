% s2let_demo7
% -----------------------------------------------------------
% Curvelets: 
% Compute and plot the harmonic tiling and the curvelet kernels.
%
% ------ About this demo -------- 
% This is the prototype MATLAB code for curvelets, 
% in which curvelets are constructed using 
% [kappa kappa0] =  s2let_transform_axisym_tiling(B, L, J_min);
% with the implementation of  
%    m = el (so to obtain the hyperbolic scaling relation);
% and the implemention of the conjugate relation: 
%   %for positive m 
%   ind_pm = ssht_elm2ind(el, m);
%   % Curvelet coefficients: 
%   flm(ind_pm) = kappa(j+1,el+1);
%   %for negative m 
%   ind_nm = ssht_elm2ind(el, -m);
%   flm(ind_nm) = (-1)^m * conj(flm(ind_pm));
% -----------------------------------------------------------
% Log: 
% -  Last modified on 6th Jan 2015 : Debug completed. 
% -----------------------------------------------------------
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details
% -----------------------------------------------------------

Spin = 0;
B = 2;   
L = 16;
% For curvelets: m=l, the azimuthal band limit N = overall band limit L
N= L; 

% Minimum and maximum scale probed by wavelets 
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


% Rotate function on the sphere such that x-axis -> z-axis  (beta =pi/2) 
% (i.e. curvelets lie on the sphere top)
alpha = pi 
beta = pi/2 
gamma = 0
% Precompute Wigner small-d functions
d = zeros(L, 2*L-1, 2*L-1);
d(1,:,:) = ssht_dl(squeeze(d(1,:,:)), L, 0, beta);
for el = 1:L-1
    d(el+1,:,:) = ssht_dl(squeeze(d(el,:,:)), L, el, beta);
end

[kappa kappa0] =  s2let_transform_axisym_tiling(B, L, J_min);
% where 
% kappa0 = scaling function with kappa0[l]= sqrt(k_lambda[1+Jmin*L])
% kappa[l+jL] = wavelet generating functions = sqrt(k_lambda[1+(j+1)*L]-k_lambda[1+j*L]) 

% Plot of tiling of curvelets in harmonic space (i.e. the same as plotting the axisymmetric wavelets):
% s2let_plot_axisym_tiling(B, L, J_min);
% fname = ['s2let_d7_', configstr, '_tiling.png']
% print('-r500', '-dpng', fname)
% movefile(fname,'/Users/jenniferyhchan/WaveletsCode_PhD/s2let/figs/tiling')


% Plot curvelets on the sphere
figure('Position',[100 100 1200 600])
flm = zeros(L^2,1);
ind_pm = 0;
ind_nm = 0;
ind=0;
for j = J_min:J
for el = 0:L-1
   m = el;
   %for positive m 
   ind_pm = ssht_elm2ind(el, m);
   % Curvelet coefficients: 
   flm(ind_pm) = kappa(j+1,el+1);
   %for negative m 
   ind_nm = ssht_elm2ind(el, -m);
   flm(ind_nm) = (-1)^m * conj(flm(ind_pm));
end

% Rotate the wavelets coefficients 
   flm_rot = ssht_rotate_flm(flm, d, alpha, gamma);

   if Spin == 0
% Compute the function (non-rotated): 
%       f = ssht_inverse(flm, L, 'Reality', true);
% check imaginary components:
% f = ssht_inverse(flm, L, 'Reality', false);
% max(abs(imag(f(:))))
% 1.1232e-14
% i.e. imaginary components are negligible 

% Compute the function (rotated): 
       f_rot = ssht_inverse(flm_rot, L, 'Reality', true);
       
       ind = ind + 1;
       if ind <= maxfigs
           h = subplot(ny, nx, ind);
% Plot the non-rotated function:            
%           ssht_plot_sphere(f, L, 'Type', type, 'Lighting', lighting);

% Plot the rotated function on the sphere
           ssht_plot_sphere(f_rot, L, 'Type', type, 'Lighting', lighting);

           title(h, ['Curvelet j = ',int2str(j-J_min+1)])
           locate = get(h,'title');
           pos = get(locate,'position');
           pos(1,2) = pos(1,2)+0.7;
           pos(1,1) = pos(1,1)-0.7;
           set(locate,'pos',pos);
           v = caxis;
           temp = max(abs(v));
           caxis([-temp temp]*plot_caxis_scale)
           zoom(zoomfactor)
% view along the x-axis 
%           view([-90,0])
       end
   end
end

colormap(jet)
%%fname = [pltroot,'/s2let_demo7_', configstr, '_cur_jet.png']
fname = ['s2let_demo7_', configstr, '_cur_jet.png']
print('-r200', '-dpng', fname)

%colormap(hot)
%%fname = [pltroot,'/s2let_demo7_', configstr, '_cur_hot.png']
%fname = ['s2let_demo7_', configstr, '_cur_hot.png']
%print('-r200', '-dpng', fname)


%% Plot scaling function which has no directional (i.e. m) dependence, 
%% the same as the axisymmetric wavelets and directional wavelets 
figure('Position',[100 100 300 300])
h=subplot(1, 1, 1);
flm = zeros(L^2,1);
for l = 0:L-1
     flm(l^2+l+1,1) = kappa0(l+1);
end
f = ssht_inverse(flm, L, 'Reality', true);

ssht_plot_sphere(f, L, 'Type', type, 'Lighting',  lighting);
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
%%fname = [pltroot,'/s2let_demo7_', configstr, '_scal_jet.png']
fname = ['s2let_demo7_', configstr, '_scal_jet.png']
print('-r200', '-dpng', fname)

% colormap(hot)
%%fname = [pltroot,'/s2let_demo7_', configstr, '_scal_hot.png']
%fname = ['s2let_demo7_', configstr, '_scal_hot.png']
%print('-r200', '-dpng', fname)
