% s2let_demo7_Curvelet_Plot_Funcs_and_Tiling
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
clear all;

% ---------------
% Define curvelet parameters: 
% ---------------
Spin = 0;  % Spin value of curvelet 
B = 2;     % B = 2 for dyadic sampling 
L = 64;    % the overall band limit
N= L;      % for curvelest, the azimuthal band limit N = L
J_min = 2; % the minimum scale probed by curvelets 
J =s2let_jmax(L, B); % the maximum scale probed by curvelets =ceil(log L/ log B);  

% ---------------
% Plot the curvelets on the sphere:
% ---------------
disp(' - Plot the curvelets on the sphere');
s2let_plot_spin_curvelet_on_sphere(B, L, J_min, ...
                                   'Spin', Spin, ...
                                   'SpinLowered', false, ...
                                   'SpinLoweredFrom', 0);


% ---------------
% Plot the tiling of scaling function and curvelets in harmonic space:
% And plot the scaling function and curvelet functions in real space:
% ---------------
disp(' - Plot the tiling of curvelets');
s2let_plot_spin_curvelet_tiling(B, L, J_min,...
                                'Spin', Spin, ...
                                'SpinLowered', false, ...
                                'SpinLoweredFrom', 0);
               


