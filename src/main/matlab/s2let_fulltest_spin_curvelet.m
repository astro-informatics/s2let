% s2let_fulltest
% Run all exactness tests for the MW sampling,
% all wavelet transforms must reconstruct the input maps
% at floating-point precision. Various parameters are tested.
%
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

clear all;
close all;

% Main parameters
L = 16;
N = L; % For curvelets: N=L  
B = 2;  %4
Spin = 0; 
J_min = 1; %1
J = s2let_jmax(L, B);

% ------------------------
% Test curvelets: 
% ------------------------
disp('Checks tiling of harmonic space for curvelets')
[cur_lm scal_l] = s2let_spin_curvelet_tiling(B, L, J_min);
error_on_cur_tiling = s2let_check_cur_tiling(cur_lm, scal_l, L, Spin, J, J_min)


% Generate signals: 
disp('-------------')
disp('Generate band-limited function')
flm_gen_cur = zeros(L^2,1);
flm_gen_cur = rand(size(flm_gen_cur)) + sqrt(-1)*rand(size(flm_gen_cur));
flm_gen_cur = 2.*(flm_gen_cur - (1+sqrt(-1))./2);
disp('Construct the corresponding signal on the sphere')
f_gen_cur = ssht_inverse(flm_gen_cur, L, 'Method', 'MW');
disp('Construct the corresponding spin signal on the sphere')
f_spin_gen_cur = ssht_inverse(flm_gen_cur, L, 'Method', 'MW', 'Spin', Spin);



disp('==============')
disp('Complex Signals, Full resolution tests start (Upsample: true):')
disp('==============')
% Test: signal reconstruction in harmonic space: 
disp('ana_lm2lmn : Perform (spin) harmonic-to-Wigner transform (lm2lmn) with custom parameters')
[f_cur_lmn, f_scal_lm]= s2let_transform_spin_curvelet_analysis_lm2lmn(flm_gen_cur, cur_lm, scal_l, ...
                                                                 'B', B, 'L', L, 'J_min', J_min,...
                                                                 'Spin', Spin,  'Reality', false, 'Upsample', true,...
                                                                 'SpinLowered', false,  'SpinLoweredFrom', 0, ...
                                                                  'Sampling', 'MW');
disp('syn_lmn2lm : Perform inverse transform (lmn2lm) with custom parameters')
flm_rec_spin_lmn2lm_custom = s2let_transform_spin_curvelet_synthesis_lmn2lm(f_cur_lmn, f_scal_lm, cur_lm, scal_l,...
                                                                      'B', B, 'L', L, 'J_min', J_min,...
                                                                      'Spin', Spin,  'Reality', false, 'Upsample', true,...
                                                                      'SpinLowered', false,  'SpinLoweredFrom', 0, ...
                                                                      'Sampling', 'MW');

default = max(abs(flm_gen_cur-flm_rec_spin_lmn2lm_custom))

disp('-------------')
% Test: signal reconstruction in pixel space: 
disp('ana_lm2cur: Perform (spin) harmonic-to-curvelet transform (lm2cur) with custom parameters')
[f_cur, f_scal] = s2let_transform_spin_curvelet_analysis_lm2cur(flm_gen_cur,  'B', B, 'L', L, 'J_min', J_min, 'Spin', Spin, 'Upsample',true);
disp('syn_cur2lm:  Perform inverse transform (cur2lm) with custom parameters')
flm_rec_spin_lm2cur_custom = s2let_transform_spin_curvelet_synthesis_cur2lm(f_cur, f_scal, 'B', B, 'L', L, 'J_min', J_min, 'Spin', Spin, 'Upsample', true);
default = max(abs(flm_gen_cur-flm_rec_spin_lm2cur_custom))

%disp('ana_lm2cur: Perform scalar harmonic-to-curvelet transform (lm2cur) with default parameters')
%[f_cur, f_scal] = s2let_transform_curvelet_analysis_lm2cur(flm_gen_cur, 'Upsample', true);
%flm_rec_lm2cur = s2let_transform_curvelet_synthesis_cur2lm(f_cur, f_scal, 'Upsample', true);
%default = max(abs(flm_gen_cur-flm_rec_lm2cur))

%disp('Perform (spin) harmonic-to-curvelet transform (lm2cur) with default parameters')
%[f_cur, f_scal] = s2let_transform_curvelet_analysis_lm2cur(flm_gen_cur, 'N', N, 'Spin', Spin, 'Upsample', true);
%flm_rec_spin_lm2cur = s2let_transform_curvelet_synthesis_cur2lm(f_cur, f_scal, 'N', N, 'Spin', Spin, 'Upsample', true);
%default = max(abs(flm_gen_cur-flm_rec_spin_lm2cur))


disp('-------------')
% Test: signal reconstruction in pixel space: 
disp('ana_px2cur: Perform (spin) curvelet transform (px2cur) with custom parameters')
[f_cur, f_scal] = s2let_transform_spin_curvelet_analysis_px2cur(f_spin_gen_cur ,...
                                                    'B', B, 'L', L, ...
                                                    'J_min', J_min, 'Spin', Spin, ...
                                                    'Upsample', true, 'Reality', false,...
                                                    'SpinLowered', false, ...
                                                    'SpinLoweredFrom', 0,...
                                                    'Sampling', 'MW');
disp('syn_cur2px:  Perform inverse transform (cur2px) with custom parameters')
f_rec_spin_px2cur_custom =  s2let_transform_spin_curvelet_synthesis_cur2px(f_cur, f_scal,  ...
                                                   'B', B, 'L', L, ...
                                                    'J_min', J_min, 'Spin', Spin, ...
                                                    'Upsample', true, 'Reality', false,...
                                                    'SpinLowered', false, ...
                                                    'SpinLoweredFrom', 0,...
                                                    'Sampling', 'MW');

default = max(abs(f_spin_gen_cur (:)-f_rec_spin_px2cur_custom(:)))

%disp('Perform (scalar) curvelet transform with default parameters')
%[f_cur, f_scal] = s2let_transform_analysis_cur_mw(f_gen_cur, 'N', N, 'Upsample', true);
%f_rec_usingcur_syn = s2let_transform_synthesis_cur_mw(f_cur, f_scal, 'N', N, 'Upsample', true);
%default = max(max(abs(f_gen_cur-f_rec_usingcur_syn)))

%disp('Perform (spin) curvelet transform with default parameters')
%[f_cur, f_scal] = s2let_transform_analysis_cur_mw(f_s_cur, 'N', N, 'Spin', Spin, 'Upsample', true);
%f_rec_usingcur = s2let_transform_synthesis_cur_mw(f_cur, f_scal, 'N', N, 'Spin', Spin, 'Upsample', true);
%default = max(max(abs(f_s_cur-f_rec_usingcur)))


disp('==============')
disp('Complex Signals, multi-resolution tests start (Upsample: false):')
disp('==============')

% Test: signal reconstruction in harmonic space: 
disp('ana_lm2lmn : Perform (spin) harmonic-to-Wigner transform (lm2lmn) with custom parameters')
[f_cur_lmn, f_scal_lm]= s2let_transform_curvelet_analysis_lm2lmn(flm_gen_cur, cur_lm, scal_l, ...
                                                                 'B', B, 'L', L, 'J_min', J_min,...
                                                                 'Spin', Spin, ...
                                                                 'Reality', false, ...
                                                                 'Upsample', false,...
                                                                 'SpinLowered', false,  'SpinLoweredFrom', 0, ...
                                                                  'Sampling', 'MW');
disp('syn_lmn2lm : Perform inverse transform (lmn2lm) with custom parameters')
flm_rec_spin_lmn2lm_custom = s2let_transform_curvelet_synthesis_lmn2lm(f_cur_lmn, f_scal_lm, cur_lm, scal_l,...
                                                                      'B', B, 'L', L, 'J_min', J_min,...
                                                                      'Spin', Spin, ...
                                                                      'Reality', false, ...
                                                                      'Upsample', false,...
                                                                      'SpinLowered', false,  'SpinLoweredFrom', 0, ...
                                                                      'Sampling', 'MW');

default = max(abs(flm_gen_cur-flm_rec_spin_lmn2lm_custom))

disp('-------------')
% Test: signal reconstruction in pixel space: 
disp('ana_lm2cur: Perform (spin) harmonic-to-curvelet transform (lm2cur) with custom parameters')
[f_cur, f_scal] = s2let_transform_curvelet_analysis_lm2cur(flm_gen_cur,  'B', B, 'L', L, 'J_min', J_min, 'Spin', Spin, 'Upsample', false);
disp('syn_cur2lm:  Perform inverse transform (cur2lm) with custom parameters')
flm_rec_spin_lm2cur_custom = s2let_transform_curvelet_synthesis_cur2lm(f_cur, f_scal, 'B', B, 'L', L, 'J_min', J_min, 'Spin', Spin, 'Upsample', false);
default = max(abs(flm_gen_cur-flm_rec_spin_lm2cur_custom))

disp('-------------')
% Test: signal reconstruction in pixel space: 
disp('ana_px2cur: Perform (spin) curvelet transform (px2cur) with custom parameters')
[f_cur, f_scal] = s2let_transform_curvelet_analysis_px2cur(f_spin_gen_cur ,...
                                                    'B', B, 'L', L, ...
                                                    'J_min', J_min, 'Spin', Spin, ...
                                                    'Upsample', false, 'Reality', false,...
                                                    'SpinLowered', false, ...
                                                    'SpinLoweredFrom', 0,...
                                                    'Sampling', 'MW');
disp('syn_cur2px:  Perform inverse transform (cur2px) with custom parameters')
f_rec_spin_px2cur_custom =  s2let_transform_curvelet_synthesis_cur2px(f_cur, f_scal,  ...
                                                   'B', B, 'L', L, ...
                                                    'J_min', J_min, 'Spin', Spin, ...
                                                    'Upsample', false, 'Reality', false,...
                                                    'SpinLowered', false, ...
                                                    'SpinLoweredFrom', 0,...
                                                    'Sampling', 'MW');

default = max(abs(f_spin_gen_cur (:)-f_rec_spin_px2cur_custom(:)))


%disp('Perform multiresolution (spin) curvelet transform with default parameters')
%[f_cur, f_scal] = s2let_transform_analysis_cur_mw(f_s_cur, 'N', N, 'Spin', Spin);
%f_rec_usingcur = s2let_transform_synthesis_cur_mw(f_cur, f_scal, 'N', N, 'Spin', Spin);
%default = max(max(abs(f_s_cur-f_rec_usingcur)))


%disp('Perform (spin) curvelet transform (analysis and synthesis) with custom parameters')
%[f_cur, f_scal] = s2let_transform_analysis_cur_mw(f_s_cur,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', false);
%f_rec_usingcur = s2let_transform_synthesis_cur_mw(f_cur, f_scal,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', false);
%default = max(max(abs(f_s_cur-f_rec_usingcur)))

%disp('Perform multiresolution (scalar) curvelet transform with default parameters')
%[f_cur, f_scal] = s2let_transform_analysis_cur_mw(f_gen_cur, 'N', N);
%f_rec_usingcur = s2let_transform_synthesis_cur_mw(f_cur, f_scal, 'N', N);
%default = max(max(abs(f_gen_cur-f_rec_usingcur)))


%disp('Perform multiresolution (spin) curvelet transform(analysis and synthesis) with custom parameters')
%[f_cur, f_scal] = s2let_transform_analysis_cur_mw(f_s_cur,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', true);
%f_rec_usingcur = s2let_transform_synthesis_cur_mw(f_cur, f_scal,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', true);
%default = max(max(abs(f_s_cur-f_rec_usingcur)))


%{
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
f_real_cur = ssht_inverse(flm_gen_cur, L, 'Method', 'MW', 'Reality', true);


disp('Perform multiresolution real curvelet transform (analysis and synthesis) with default parameters')
[f_cur, f_scal] = s2let_transform_analysis_cur_mw(f_real_cur, 'N', N, 'Reality', true);
f_rec_usingcur = s2let_transform_synthesis_cur_mw(f_cur, f_scal, 'N', N, 'Reality', true);
default = max(max(abs(f_real_cur-f_rec_usingcur)))

%disp('Perform full resolution real curvelet transform with default parameters')
%[f_cur, f_scal] = s2let_transform_analysis_cur_mw(f_real_cur, 'N', N, 'Reality', true, 'Upsample', true);
%f_rec_usingcur = s2let_transform_synthesis_cur_mw(f_cur, f_scal, 'N', N, 'Reality', true, 'Upsample', true);
%default = max(max(abs(f_real_cur-f_rec_usingcur)))

%}



