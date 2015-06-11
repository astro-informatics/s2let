% s2let_demo8
% ------ About this demo -------- 
% 1a) Compute and plot the harmonic tiling and the wavelet kernels.
% 1b) Plot the scaling function and curvelet functions in real space. 
% 2)  Plot the curvelets on the Sphere. 
% 
% Curvelets and the scaling function are generated via the matlab function 
%  "s2let_curvelet_tiling(B, L, J_min, <options>)".  
% External matlab function  
%  "s2let_plot_cur_tiling(B, L, J_min, 'spin', Spin);" is used to accomplish 1a) and 1b) 
%  "s2let_plot_cur_on_sphere" is called to accomplish 2). 
% -----------------------------------------------------------
% Log: 
% -  constructed by Jennifer Y H Chan on 5th June 2014  
% -----------------------------------------------------------
% S2LET package to perform wavelet transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details
% 
% Modified S2LET package to perform Curvelets on the Sphere.
% -----------------------------------------------------------


% ---------------
% Define curvelet parameters: 
% ---------------
Spin = 0;
B = 2;   % for dyadic sampling 
L = 64;
N= L;     % Since m=l, the azimuthal band limit N = overall band limit L
J_min = 2; % minimum and maximum scale probed by wavelets 
J =s2let_jmax(L, B); %=ceil(log L/ log B);  

% ---------------
% Plot the tiling of scaling function and curvelets in harmonic space:
% And plot the scaling function and curvelet functions in real space:
% ---------------
disp(' - Plot the tiling of curvelets');
s2let_plot_cur_tiling(B, L, J_min,...
                     'Spin', Spin, ...
                     'SpinLowered', false, ...
                     'SpinLoweredFrom', 0);

% ---------------
% Plot the curvelets on the sphere:
% Define Euler angles (for rotation):
% ---------------
alpha =  pi ;
beta = pi/2 ;
gamma = 0 ;
disp(' - Plot the curvelets on the sphere');
s2let_plot_cur_on_sphere(alpha, beta, gamma, B, L, J_min, ...
                         'Spin', Spin, ...
                         'SpinLowered', false, ...
                         'SpinLoweredFrom', 0);


