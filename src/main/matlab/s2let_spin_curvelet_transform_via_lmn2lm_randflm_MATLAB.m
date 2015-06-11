% s2let_curvelet_matlab_prototype
% Run curvelet analysis and synthesis 
% of randomly generated signals f 
%
%
% Options consist of parameter type and value pairs.
% Valid options include:
%
%  'Upsample'      = { false   [multiresolution algorithm (default)],
%                      true    [full resolution wavelets] },
% -----------------------------------------------------------
% Log: 
% -  constructed by Jennifer Y H Chan on 5th June 2015  
% -----------------------------------------------------------
% S2LET package to perform wavelet transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details
% -----------------------------------------------------------

clear all;
close all;

% Main parameters
Spin = 0;
L = 16;
N = L; % For curvelets: N=L  %N=1 for axisymmetric curvelets
B = 2;
J_min = 2;  
J =s2let_jmax(L, B);  %=ceil(log L/ log B);  

%{
disp('Generates random band-limited function')
flm_gen = zeros(L^2,1);
flm_gen = rand(size(flm_gen)) + sqrt(-1)*rand(size(flm_gen));
flm_gen = 2.*(flm_gen - (1+sqrt(-1))./2);

disp('Construct the corresponding signal on the sphere')
f_gen = ssht_inverse(flm_gen, L, 'Method', 'MW');
%}

flm_gen = zeros(L^2,1);
disp('read from file flm_gen'); 
fid= fopen('/Users/jenniferyhchan/WaveletsCode_PhD/s2let_curvelets_MATLAB/1_cur_flm_randgen_mw_test.dat');
rawData=fscanf(fid, '%f, %f',[2 256]);
fclose(fid);
complexData=complex(rawData(1,:),rawData(2,:));
flm_gen= complexData.' ;  % Non-conjugate transpose
disp('Construct the corresponding signal via ssht_inverse'); 
f_gen = ssht_inverse(flm_gen, L, 'Method', 'MW');
flm_gen= ssht_forward(f_gen, L, 'Method', 'MW');

% ---------------
% Tile curvelets:
% ---------------
% Call curvelet- and scaling-function- generating functions 
disp('curvelet_tiling: Tile curvelets in harmonic space (cur_lm, scal_l)')
[cur_lm scal_l] = s2let_spin_curvelet_tiling(B, L, J_min, ...
                                            'Spin', Spin, 'SpinLowered', false,  'SpinLoweredFrom', 0);
% Check tiling error:
error_on_cur_tiling = s2let_check_cur_tiling(cur_lm, scal_l, L, Spin, J, J_min)

% -----------------
% Signal analysis: (harmonic to wigner space) 
% -----------------
% Call matlab function 's2let_transform_curvelet_analysis_lm2lmn'
disp('analysis_lm2lmn...')
[f_cur_lmn, f_scal_lm] = s2let_transform_spin_curvelet_analysis_lm2lmn(flm_gen, cur_lm, scal_l, ...
                                                         'B', B, 'L', L, 'J_min', J_min,...
                                                         'Spin', Spin, ...
                                                         'Reality', false, ...
                                                         'Upsample', false,...
                                                         'SpinLowered', false,  'SpinLoweredFrom', 0, ...
                                                         'Sampling', 'MW');
% -----------------
% Signal synthesis: (pixel to harmonic space) 
% -----------------
% Call matlab function 's2let_transform_curvelet_synthesis_lmn2lm'
disp('synthesis_lmn2lm...')
flm_rec  = s2let_transform_spin_curvelet_synthesis_lmn2lm(f_cur_lmn, f_scal_lm, cur_lm, scal_l, ...
                                                      'B', B, 'L', L, 'J_min', J_min,...
                                                      'Spin', Spin, ...
                                                      'Reality', false,...
                                                      'Upsample', false,...
                                                      'SpinLowered', false,  'SpinLoweredFrom', 0, ...
                                                      'Sampling', 'MW');
                                       
disp('Compute the re-constructed function via ssht_inverse ');
f_rec = ssht_inverse(flm_rec, L, 'Method', 'MW'); 

disp('Both analysis and synthesis are done! ');
disp('');
disp('- Test exact transform: check the difference between flm_gen and flm_rec:');
maxerr = max(abs(flm_gen - flm_rec))
disp('Check the difference between f_gen and f_rec: ');
maxerr = max(abs(f_gen(:) - f_rec(:)))