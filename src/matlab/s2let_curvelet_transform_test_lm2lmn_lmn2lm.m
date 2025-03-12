% s2let_curvelet_transform_test_lm2lmn_lmn2lm
%
% Run curvelet analysis (harmonic to Wigner space) and 
% synthesis (Wigner to harmonic space)
% of randomly generated signals f and check exactness
%
% Options consist of parameter type and value pairs.
% Valid options include:
%  'B'               = { Dilation factor; B > 1 (default = 2) }
%  'L'               = { Harmonic band-limit; L > 0 (default = Lguessed) }
%  'J_min'           = { the minimal wavelet scale,(default = 0)}
%  'Spin'            = { Spin number; Spin >= 0 (default = 0) }
%  'Reality'         = { false   [do not assume corresponding signal f real (default)],
%                        true    [assume f real (improves performance)] }
%  'Upsample'        = { false   [multiresolution algorithm (default)],
%                        true    [full resolution wavelets] },
% -----------------------------------------------------------
% S2LET package to perform Wavelet Transform on the Sphere.
% Copyright (C) 2012-2016  Boris Leistedt, Jennifer Chan & Jason McEwen
% See LICENSE.txt for license details
% -----------------------------------------------------------

clear all ;
close all;

% Curvelet parameters
Spin = 0;             % Spin value of curvelet 
L = 16;               % Angular band-limit
B = 2;                % B = 2 for dyadic sampling
J_min = 2;            % Minimum scale probed by curvelets
J =s2let_jmax(L, B);  % Maximum scale probed by curvelets =ceil(log L/ log B);


disp('Generates random band-limited function')
flm_gen = zeros(L^2,1);
flm_gen = rand(size(flm_gen)) + sqrt(-1)*rand(size(flm_gen));
flm_gen = 2.*(flm_gen - (1+sqrt(-1))./2);
disp('Construct the corresponding signal on the sphere')
f_gen = ssht_inverse(flm_gen, L,'Method', 'MW');
disp('Construct the corresponding spin signal on the sphere')
f_spin_gen = ssht_inverse(flm_gen, L, 'Spin', Spin, 'Method', 'MW');
flm_spin_gen= ssht_forward(f_spin_gen, L, 'Spin', Spin, 'Method', 'MW');
disp('----------- ');


% ---------------
% Tile curvelets:
% ---------------
disp('curvelet_tiling: Tile curvelets in harmonic space (cur_lm, scal_l)')
% Call curvelet- and scaling-function- generating functions 
[cur_lm scal_l] = s2let_curvelet_tiling(B, L, J_min, ...
                                            'Spin', Spin, 'SpinLowered', false,  'SpinLoweredFrom', 0);
% Check tiling error: 
disp('Check if admissibility condition is satisfied: ');  
error_on_cur_tiling = s2let_check_cur_tiling(cur_lm, scal_l, L, Spin, J, J_min)
disp('----------- ');  


disp(' ');
% ================== FULL-RESOLUTION ===================% 
disp('Curvelet transform: Full resolution (Upsample: true)');
% -----------------
% Signal analysis: (harmonic to Wigner space) 
% -----------------
disp('Spin signal, Full resolution: analysis_lm2lmn...')
[f_cur_lmn, f_scal_lm] = s2let_transform_curvelet_analysis_lm2lmn(flm_spin_gen, cur_lm, scal_l, ...
                                                                  'B', B, 'L', L, 'J_min', J_min,...
                                                                  'Spin', Spin, ...
                                                                  'Reality', false, ...
                                                                  'Upsample', true,...
                                                                  'SpinLowered', false, ...
                                                                  'SpinLoweredFrom', 0);
% -----------------
% Signal synthesis: (Wigner to harmonic space) 
% -----------------
disp('Spin signal, Full resolution: synthesis_lmn2lm...')
flm_spin_rec  = s2let_transform_curvelet_synthesis_lmn2lm(f_cur_lmn, f_scal_lm, cur_lm, scal_l, ...
                                                         'B', B, 'L', L, 'J_min', J_min,...
                                                         'Spin', Spin, ...
                                                         'Reality', false,...
                                                         'Upsample', true,...
                                                         'SpinLowered', false, ...
                                                         'SpinLoweredFrom', 0);
                                       
disp('Compute the re-constructed function via ssht_inverse ');
f_spin_rec = ssht_inverse(flm_spin_rec, L, 'Spin', Spin, 'Method', 'MW');
disp('- Test exact transform:');
disp('Check the difference between flm_gen and flm_rec:');
maxerr = max(abs(flm_spin_gen - flm_spin_rec))
disp('Check the difference between f_gen and f_rec: ');
maxerr = max(abs(f_spin_gen(:) - f_spin_rec(:)))
disp('----------- ');           

% ================== MULTI-RESOLUTION ===================% 
disp('Curvelet transform: multiresolution (Upsample: false): ');
% -----------------
% Signal analysis: (harmonic to wigner space) 
% -----------------
disp('Spin signal, Multi-resolution: analysis_lm2lmn...')
[f_cur_lmn, f_scal_lm] = s2let_transform_curvelet_analysis_lm2lmn(flm_spin_gen, cur_lm, scal_l, ...
                                                                  'B', B, 'L', L, 'J_min', J_min,...
                                                                  'Spin', Spin, ...
                                                                  'Reality', false, ...
                                                                  'Upsample', false,...
                                                                  'SpinLowered', false, ...
                                                                  'SpinLoweredFrom', 0);
% -----------------
% Signal synthesis: (pixel to harmonic space) 
% -----------------
disp('Spin signal, Multi-resolution: synthesis_lmn2lm...')
flm_spin_rec  = s2let_transform_curvelet_synthesis_lmn2lm(f_cur_lmn, f_scal_lm, cur_lm, scal_l, ...
                                                          'B', B, 'L', L, 'J_min', J_min,...
                                                          'Spin', Spin, ...
                                                          'Reality', false,...
                                                          'Upsample', false,...
                                                          'SpinLowered', false, ...
                                                          'SpinLoweredFrom', 0);

disp('Compute the re-constructed function via ssht_inverse ');
f_spin_rec = ssht_inverse(flm_spin_rec, L, 'Spin', Spin, 'Method', 'MW');
disp('- Test exact transform:');
disp('Check the difference between flm_gen and flm_rec:');
maxerr = max(abs(flm_spin_gen - flm_spin_rec))
disp('Check the difference between f_gen and f_rec: ');
maxerr = max(abs(f_spin_gen(:) - f_spin_rec(:)))
disp('----------- ');


disp(' ')
disp('=============================')
disp('REAL Signals TEST')
disp('=============================')
% -------------------------
disp('Constraint on flms to generate real signal')
% -------------------------
for el = 0:L-1
    ind = el*el + el + 1;
    flm_gen(ind,1) = real(flm_gen(ind,1));
    for m = 1:el
       ind_pm = el*el + el + m + 1;
       ind_nm = el*el + el - m + 1;
       flm_gen(ind_nm,1) = (-1)^m * conj(flm_gen(ind_pm,1));
    end
end

% ---------------
% Tile spin-0 curvelets:
% ---------------
disp('Real signal: curvelet_tiling: Tile spin-0 curvelets in harmonic space (cur_lm, scal_l)')
[cur_lm scal_l] = s2let_curvelet_tiling(B, L, J_min, 'Spin', 0);
disp('----------- ');    

% ================== FULL-RESOLUTION ===================% 
disp('Curvelet transform: Full-resolution (Upsample: true): ');
% -----------------
% Signal analysis: (harmonic to Wigner space) 
% -----------------
disp('Real signal, Full resolution: analysis_lm2lmn...')
[f_cur_lmn, f_scal_lm] = s2let_transform_curvelet_analysis_lm2lmn(flm_gen, cur_lm, scal_l, ...
                                                                  'B', B, 'L', L, 'J_min', J_min,...
                                                                  'Reality', true, ...
                                                                  'Upsample', true);
% -----------------
% Signal synthesis: (Wigner to harmonic space) 
% -----------------
disp('Real signal, Full resolution: synthesis_lmn2lm...')
flm_rec  = s2let_transform_curvelet_synthesis_lmn2lm(f_cur_lmn, f_scal_lm, cur_lm, scal_l, ...
                                                     'B', B, 'L', L, 'J_min', J_min,...
                                                     'Reality', true,...
                                                     'Upsample', true);
f_gen = ssht_inverse(flm_gen, L,'Method', 'MW');
f_rec = ssht_inverse(flm_rec, L,'Method', 'MW');
disp('- Test exact transform:');
disp('Check the difference between flm_gen and flm_rec:');
maxerr = max(abs(flm_gen - flm_rec))
disp('Check the difference between f_gen and f_rec: ');
maxerr = max(abs(f_gen(:) - f_rec(:)))
disp('----------- ');          

% ================== MULTI-RESOLUTION ===================% 
disp('Curvelet transform: Multi-resolution (Upsample: false): ');
% -----------------
% Signal analysis: (harmonic to Wigner space) 
% -----------------
disp('Real signal, Multi-reolution: analysis_lm2lmn...')
[f_cur_lmn, f_scal_lm] = s2let_transform_curvelet_analysis_lm2lmn(flm_gen, cur_lm, scal_l, ...
                                                                  'B', B, 'L', L, 'J_min', J_min,...
                                                                  'Reality', true, ...
                                                                  'Upsample', false);
% -----------------
% Signal synthesis: (Wigner to harmonic space) 
% -----------------
disp('Real signal, Multi-reolution: synthesis_lmn2lm...')
flm_rec  = s2let_transform_curvelet_synthesis_lmn2lm(f_cur_lmn, f_scal_lm, cur_lm, scal_l, ...
                                                     'B', B, 'L', L, 'J_min', J_min,...
                                                     'Reality', true,...
                                                     'Upsample', false);
f_gen = ssht_inverse(flm_gen, L,'Method', 'MW');
f_rec = ssht_inverse(flm_rec, L,'Method', 'MW');
disp('- Test exact transform:');
disp('Check the difference between flm_gen and flm_rec:');
maxerr = max(abs(flm_gen - flm_rec))
disp('Check the difference between f_gen and f_rec: ');
maxerr = max(abs(f_gen(:) - f_rec(:)))
disp('----------- ');