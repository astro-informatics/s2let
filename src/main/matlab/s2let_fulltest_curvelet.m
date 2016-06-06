% s2let_fulltest
% Run all exactness tests for the MW sampling,
% all wavelet transforms must reconstruct the input maps
% at floating-point precision. Various parameters are tested.
%
% -----------------------------------------------------------
% S2LET package to perform Wavelet Transform on the Sphere.
% Copyright (C) 2012-2016  Boris Leistedt, Jennifer Chan & Jason McEwen
% See LICENSE.txt for license details
% -----------------------------------------------------------

clear all;
close all;

% Curvelet parameters
L = 64;
B = 2; 
Spin = 0;  
J_min = 1; 
J = s2let_jmax(L, B);
disp('-------------') 
% ------------------------
% Tile curvelets: 
% ------------------------
disp('Checks the harmonic tiling of curvelets')
[cur_lm scal_l] = s2let_curvelet_tiling(B, L, J_min, ...
                                        'Spin', Spin, 'SpinLowered', false,  'SpinLoweredFrom', 0);
error_on_cur_tiling = s2let_check_cur_tiling(cur_lm, scal_l, L, Spin, J, J_min)


% -------------------------
disp('Generate random band-limited (complex) function')
% -------------------------
flm_gen_cur = zeros(L^2,1);
flm_gen_cur = rand(size(flm_gen_cur)) + sqrt(-1)*rand(size(flm_gen_cur));
flm_gen_cur = 2.*(flm_gen_cur - (1+sqrt(-1))./2);
disp('Construct the corresponding signal on the sphere')
f_gen_cur = ssht_inverse(flm_gen_cur, L, 'Method', 'MW');
% -------------------------
disp('Construct the corresponding spin signal on the sphere')
% -------------------------
f_spin_gen_cur = ssht_inverse(flm_gen_cur, L, 'Spin', Spin, 'Method', 'MW');
flm_spin_gen_cur = ssht_forward(f_spin_gen_cur, L, 'Spin', Spin, 'Method', 'MW');


disp('==============')
disp('Complex Signals, Full resolution tests start (Upsample: true):')
disp('==============')
disp('ana_lm2lmn : Perform (spin) harmonic-to-Wigner transform (lm2lmn) with custom parameters')
[f_cur_lmn, f_scal_lm]= s2let_transform_curvelet_analysis_lm2lmn(flm_spin_gen_cur, cur_lm, scal_l, ...
                                                                 'B', B, 'L', L, 'J_min', J_min,...
                                                                 'Spin', Spin, ...
                                                                 'Reality', false,...
                                                                 'Upsample', true,...
                                                                 'SpinLowered', false, ...
                                                                 'SpinLoweredFrom', 0);
disp('syn_lmn2lm : Perform inverse transform (lmn2lm) with custom parameters')
flm_rec_spin_lmn2lm_custom = s2let_transform_curvelet_synthesis_lmn2lm(f_cur_lmn, f_scal_lm, cur_lm, scal_l,...
                                                                      'B', B, 'L', L, 'J_min', J_min,...
                                                                      'Spin', Spin,  ...
                                                                      'Reality', false,...
                                                                      'Upsample', true,...
                                                                      'SpinLowered', false, ...
                                                                      'SpinLoweredFrom', 0);
default = max(abs(flm_spin_gen_cur - flm_rec_spin_lmn2lm_custom))

disp('-------------')
disp('ana_lm2cur: Perform (spin) harmonic-to-curvelet transform (lm2cur) with custom parameters')
[f_cur, f_scal] = s2let_transform_curvelet_analysis_lm2cur(flm_spin_gen_cur, ...
                                                           'B', B, 'L', L, 'J_min', J_min,...
                                                           'Spin', Spin, 'Upsample',true);
disp('syn_cur2lm:  Perform inverse transform (cur2lm) with custom parameters')
flm_rec_spin_lm2cur_custom = s2let_transform_curvelet_synthesis_cur2lm(f_cur, f_scal,...
                                                                       'B', B, 'L', L, 'J_min', J_min, ...
                                                                       'Spin', Spin, 'Upsample', true);
default = max(abs(flm_spin_gen_cur - flm_rec_spin_lm2cur_custom))

disp('-------------')
disp('ana_px2cur: Perform (spin) curvelet transform (px2cur) with custom parameters')
% c.f. s2let_transform_analysis_mw.m for directional wavelet
[f_cur, f_scal] = s2let_transform_curvelet_analysis_px2cur(f_spin_gen_cur ,...
                                                    'B', B, 'L', L, ...
                                                    'J_min', J_min, 'Spin', Spin, ...
                                                    'Upsample', true, 'Reality', false,...
                                                    'SpinLowered', false, ...
                                                    'SpinLoweredFrom', 0);
disp('syn_cur2px: Perform inverse transform (cur2px) with custom parameters')
% c.f. s2let_transform_synthesis_mw.m for directional wavelet
f_rec_spin_px2cur_custom =  s2let_transform_curvelet_synthesis_cur2px(f_cur, f_scal,  ...
                                                    'B', B, 'L', L, ...
                                                    'J_min', J_min, 'Spin', Spin, ...
                                                    'Upsample', true, 'Reality', false,...
                                                    'SpinLowered', false, ...
                                                    'SpinLoweredFrom', 0);
default = max(abs(f_spin_gen_cur(:)-f_rec_spin_px2cur_custom(:)))

disp('')
disp(' - CHECK also: the Spin Curvelet Transform of spin signals in MW sampling with default parapmeters:')
disp('ana_px2cur: Perform (spin) curvelet transform with default parameters')
[f_cur, f_scal] = s2let_transform_curvelet_analysis_px2cur(f_spin_gen_cur,...
                                                          'Spin', Spin, ...
                                                          'Upsample', true);
disp('syn_cur2px: Perform inverse transform (cur2px) with default parameters')
f_rec_spin_px2cur_custom =  s2let_transform_curvelet_synthesis_cur2px(f_cur, f_scal, ...
                                                                      'Spin', Spin, ...
                                                                      'Upsample', true);
default = max(abs(f_spin_gen_cur(:)-f_rec_spin_px2cur_custom(:)))

disp('')
disp(' -CHECK also: the Spin Curvelet Transform of scalar signalsin MW sampling with default parapmeters:')
disp('ana_px2cur: Perform (scalar) curvelet transform with default parameters')
[f_cur, f_scal] = s2let_transform_curvelet_analysis_px2cur(f_gen_cur,...
                                                           'Upsample', true);
disp('syn_cur2px: Perform inverse transform (cur2px) with default parameters')
f_rec_spin_px2cur =  s2let_transform_curvelet_synthesis_cur2px(f_cur, f_scal,  ...
                                                               'Upsample', true);
default = max(abs(f_gen_cur(:)-f_rec_spin_px2cur(:)))


disp('')
disp('==============')
disp('Complex Signals, Multi-resolution tests start (Upsample: false):')
disp('==============')
disp('ana_lm2lmn : Perform (spin) harmonic-to-Wigner transform (lm2lmn) with custom parameters')
[f_cur_lmn, f_scal_lm]= s2let_transform_curvelet_analysis_lm2lmn(flm_spin_gen_cur, cur_lm, scal_l, ...
                                                                 'B', B, 'L', L, 'J_min', J_min,...
                                                                 'Spin', Spin, ...
                                                                 'Reality', false,...
                                                                 'Upsample', false,...
                                                                 'SpinLowered', false, ...
                                                                 'SpinLoweredFrom', 0);
disp('syn_lmn2lm : Perform inverse transform (lmn2lm) with custom parameters')
flm_rec_spin_lmn2lm_custom = s2let_transform_curvelet_synthesis_lmn2lm(f_cur_lmn, f_scal_lm, cur_lm, scal_l,...
                                                                      'B', B, 'L', L, 'J_min', J_min,...
                                                                      'Spin', Spin,  ...
                                                                      'Reality', false,...
                                                                      'Upsample', false,...
                                                                      'SpinLowered', false, ...
                                                                      'SpinLoweredFrom', 0);
default = max(abs(flm_spin_gen_cur - flm_rec_spin_lmn2lm_custom))

disp('-------------')
disp('ana_lm2cur: Perform (spin) harmonic-to-curvelet transform (lm2cur) with custom parameters')
[f_cur, f_scal] = s2let_transform_curvelet_analysis_lm2cur(flm_spin_gen_cur, ...
                                                           'B', B, 'L', L, 'J_min', J_min,...
                                                           'Spin', Spin, 'Upsample',false);
disp('syn_cur2lm:  Perform inverse transform (cur2lm) with custom parameters')
flm_rec_spin_lm2cur_custom = s2let_transform_curvelet_synthesis_cur2lm(f_cur, f_scal,...
                                                                       'B', B, 'L', L, 'J_min', J_min, ...
                                                                       'Spin', Spin, 'Upsample', false);
default = max(abs(flm_spin_gen_cur - flm_rec_spin_lm2cur_custom))

disp('-------------')
disp('ana_px2cur: Perform (spin) curvelet transform (px2cur) with custom parameters')
[f_cur, f_scal] = s2let_transform_curvelet_analysis_px2cur(f_spin_gen_cur ,...
                                                    'B', B, 'L', L, ...
                                                    'J_min', J_min, 'Spin', Spin, ...
                                                    'Upsample', false, 'Reality', false,...
                                                    'SpinLowered', false, ...
                                                    'SpinLoweredFrom', 0);
disp('syn_cur2px: Perform inverse transform (cur2px) with custom parameters')
f_rec_spin_px2cur_custom =  s2let_transform_curvelet_synthesis_cur2px(f_cur, f_scal,  ...
                                                    'B', B, 'L', L, ...
                                                    'J_min', J_min, 'Spin', Spin, ...
                                                    'Upsample', false, 'Reality', false,...
                                                    'SpinLowered', false, ...
                                                    'SpinLoweredFrom', 0);
default = max(abs(f_spin_gen_cur(:)-f_rec_spin_px2cur_custom(:)))

disp('')
disp('- CHECK also: the Spin Curvelet Transform of spin signals in MW sampling with default parameters:')
disp('ana_px2cur: Perform (spin) curvelet transform with default parameters')
[f_cur, f_scal] = s2let_transform_curvelet_analysis_px2cur(f_spin_gen_cur,...
                                                           'Spin', Spin);
disp('syn_cur2px: Perform inverse transform (cur2px) with default parameters')
f_rec_spin_px2cur=  s2let_transform_curvelet_synthesis_cur2px(f_cur, f_scal, ...
                                                              'Spin', Spin);
default = max(abs(f_spin_gen_cur(:)-f_rec_spin_px2cur(:)))

disp('')
disp('- CHECK also: the Spin Curvelet Transform of scalar signals in MW sampling with default parameters:')
disp('ana_px2cur: Perform (scalar) curvelet transform with default parameters')
[f_cur, f_scal] = s2let_transform_curvelet_analysis_px2cur(f_gen_cur);
disp('syn_cur2px: Perform inverse transform (cur2px) with default parameters')
f_rec_spin_px2cur =  s2let_transform_curvelet_synthesis_cur2px(f_cur, f_scal);
default = max(abs(f_gen_cur(:)-f_rec_spin_px2cur(:)))


disp(' ')
disp('==============')
disp('REAL Signals TEST')
disp('==============')
% -------------------------
disp('Constraint on flms to generate real signal')
% -------------------------
for el = 0:L-1
ind = el*el + el + 1;
flm_gen_cur(ind,1) = real(flm_gen_cur(ind,1));
for m = 1:el
ind_pm = el*el + el + m + 1;
ind_nm = el*el + el - m + 1;
flm_gen_cur(ind_nm,1) = (-1)^m * conj(flm_gen_cur(ind_pm,1));
end
end
disp('Construct the corresponding real signal on the sphere')
f_real_gen = ssht_inverse(flm_gen_cur, L, 'Method', 'MW', 'Reality', true);

disp('ana_px2cur: Perform (REAL) curvelet transform (px2cur) with custom parameters')
[f_cur, f_scal] = s2let_transform_curvelet_analysis_px2cur(f_real_gen ,...
                                                    'B', B, 'L', L, ...
                                                    'J_min', J_min,  ...
                                                    'Upsample', true, ...
                                                    'Reality', true);
disp('syn_cur2px: Perform inverse transform (cur2px) with custom parameters')
f_rec_real_px2cur_custom =  s2let_transform_curvelet_synthesis_cur2px(f_cur, f_scal,  ...
                                                    'B', B, 'L', L, ...
                                                    'J_min', J_min, ...
                                                    'Upsample', true, ...
                                                    'Reality',true);
default = max(abs(f_real_gen(:)- f_rec_real_px2cur_custom(:)))



disp('')
disp('==============')
disp('Rea, Signals, Multi-resolution tests start (Upsample: false):')
disp('==============')
disp('ana_px2cur: Perform (REAL) curvelet transform (px2cur) with custom parameters')
[f_cur, f_scal] = s2let_transform_curvelet_analysis_px2cur(f_real_gen ,...
                                                    'B', B, 'L', L, ...
                                                    'J_min', J_min,  ...
                                                    'Upsample', false, ...
                                                    'Reality', true);
disp('syn_cur2px: Perform inverse transform (cur2px) with custom parameters')
f_rec_real_px2cur_custom =  s2let_transform_curvelet_synthesis_cur2px(f_cur, f_scal,  ...
                                                    'B', B, 'L', L, ...
                                                    'J_min', J_min, ...
                                                    'Upsample', false, ...
                                                    'Reality',true);
default = max(abs(f_real_gen(:)- f_rec_real_px2cur_custom(:)))


disp('')
disp('==============')
disp('Real Signals, DEFAULT parameters:')
disp('==============')
disp('- CHECK also: the REAL Curvelet Transform of scalar signals in MW sampling with default parameters:')
disp('ana_px2cur: Perform (scalar) curvelet transform with default parameters')
[f_cur, f_scal] = s2let_transform_curvelet_analysis_px2cur(f_real_gen ,'Reality', true);
disp('syn_cur2px: Perform inverse transform (cur2px) with default parameters')
f_rec_real_px2cur_custom =  s2let_transform_curvelet_synthesis_cur2px(f_cur, f_scal, 'Reality',true);
default = max(abs(f_real_gen(:)- f_rec_real_px2cur_custom(:)))




