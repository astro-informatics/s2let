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
J_min = 0; %1
J = s2let_jmax(L, B);

% ------------------------
% Test curvelets: 
% ------------------------
disp('Checks tiling of harmonic space for axysimmetric curvelets')
[kappa kappa0] = s2let_transform_axisym_tiling(B, L, J_min);
error_on_axisym_tiling = s2let_check_axisym_tiling(kappa, kappa0, L, J)


%disp('Checks tiling of slm for directional wavelets')
%[slm] = s2let_wavelet_tiling(B, L, N, Spin, J_min);
%error_on_direction_tiling = s2let_check_slm(slm,L)


% Curvelets:
disp('Checks tiling of harmonic space for curvelets')
[cur_lm scal_l] = s2let_curvelet_tiling(B, L, J_min);
error_on_cur_tiling = s2let_check_cur_tiling(scal_l, cur_lm, L, Spin, J, J_min)

disp('Generates band-limited function')
flm_gen_cur = zeros(L^2,1);
flm_gen_cur = rand(size(flm_gen_cur)) + sqrt(-1)*rand(size(flm_gen_cur));
flm_gen_cur = 2.*(flm_gen_cur - (1+sqrt(-1))./2);
disp('Construct the corresponding signal on the sphere')
f_gen_cur = ssht_inverse(flm_gen_cur, L, 'Method', 'MW');
disp('Construct the corresponding spin signal on the sphere')
f_s_cur = ssht_inverse(flm_gen_cur, L, 'Method', 'MW', 'Spin', Spin);


%disp('Perform scalar harmonic-to-curvelet transform (lm2cur) with default parameters')
%[f_cur, f_scal] = s2let_transform_analysis_lm2cur(flm_gen_cur, 'N', N, 'Upsample', true);
%flm_rec_lm2cur = s2let_transform_synthesis_lm2cur(f_cur, f_scal, 'N', N, 'Upsample', true);
%default = max(abs(flm_gen_cur-flm_rec_lm2cur))

%disp('Perform (spin) harmonic-to-curvelet transform (lm2cur) with default parameters')
%[f_cur, f_scal] = s2let_transform_analysis_lm2cur(flm_gen_cur, 'N', N, 'Spin', Spin, 'Upsample', true);
%flm_rec_spin_lm2cur = s2let_transform_synthesis_lm2cur(f_cur, f_scal, 'N', N, 'Spin', Spin, 'Upsample', true);
%default = max(abs(flm_gen_cur-flm_rec_spin_lm2cur))

disp('Perform (spin) harmonic-to-curvelet transform (lm2cur) with custom parameters')
[f_cur, f_scal] = s2let_transform_analysis_lm2cur(flm_gen_cur,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', false);
flm_rec_spin_lm2cur_custom = s2let_transform_synthesis_lm2cur(f_cur, f_scal,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', false);
default = max(abs(flm_gen_cur-flm_rec_spin_lm2cur_custom))

%disp('Perform (scalar) curvelet transform with default parameters')
%[f_cur, f_scal] = s2let_transform_analysis_cur_mw(f_gen_cur, 'N', N, 'Upsample', true);
%f_rec_usingcur_syn = s2let_transform_synthesis_cur_mw(f_cur, f_scal, 'N', N, 'Upsample', true);
%default = max(max(abs(f_gen_cur-f_rec_usingcur_syn)))

%disp('Perform (spin) curvelet transform with default parameters')
%[f_cur, f_scal] = s2let_transform_analysis_cur_mw(f_s_cur, 'N', N, 'Spin', Spin, 'Upsample', true);
%f_rec_usingcur = s2let_transform_synthesis_cur_mw(f_cur, f_scal, 'N', N, 'Spin', Spin, 'Upsample', true);
%default = max(max(abs(f_s_cur-f_rec_usingcur)))

disp('Perform (spin) curvelet transform (analysis and synthesis) with custom parameters')
[f_cur, f_scal] = s2let_transform_analysis_cur_mw(f_s_cur,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', false);
f_rec_usingcur = s2let_transform_synthesis_cur_mw(f_cur, f_scal,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', false);
default = max(max(abs(f_s_cur-f_rec_usingcur)))

%disp('Perform multiresolution (scalar) curvelet transform with default parameters')
%[f_cur, f_scal] = s2let_transform_analysis_cur_mw(f_gen_cur, 'N', N);
%f_rec_usingcur = s2let_transform_synthesis_cur_mw(f_cur, f_scal, 'N', N);
%default = max(max(abs(f_gen_cur-f_rec_usingcur)))

%disp('Perform multiresolution (spin) curvelet transform with default parameters')
%[f_cur, f_scal] = s2let_transform_analysis_cur_mw(f_s_cur, 'N', N, 'Spin', Spin);
%f_rec_usingcur = s2let_transform_synthesis_cur_mw(f_cur, f_scal, 'N', N, 'Spin', Spin);
%default = max(max(abs(f_s_cur-f_rec_usingcur)))

disp('Perform multiresolution (spin) curvelet transform(analysis and synthesis) with custom parameters')
[f_cur, f_scal] = s2let_transform_analysis_cur_mw(f_s_cur,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', true);
f_rec_usingcur = s2let_transform_synthesis_cur_mw(f_cur, f_scal,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', true);
default = max(max(abs(f_s_cur-f_rec_usingcur)))



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

stop


% ---------------------
% Directional wavelets:
% ---------------------
disp('')
disp('Directional wavelets')
%disp('Checks tiling of harmonic space for axysimmetric wavelets')
%[kappa kappa0] = s2let_transform_axisym_tiling(B, L, J_min);
%error_on_axisym_tiling = s2let_check_axisym_tiling(kappa, kappa0, L, J)

disp('Checks tiling of harmonic space for directional wavelets')
[psi phi] = s2let_wavelet_tiling(B, L, N, Spin, J_min);
error_on_tiling = s2let_check_tiling(psi, phi, L, Spin, J)

disp('Generates band-limited function')
flm = zeros(L^2,1);
flm = rand(size(flm)) + sqrt(-1)*rand(size(flm));
flm = 2.*(flm - (1+sqrt(-1))./2);
disp('Construct the corresponding signal on the sphere')
f = ssht_inverse(flm, L, 'Method', 'MW');
disp('Construct the corresponding spin signal on the sphere')
f_s = ssht_inverse(flm, L, 'Method', 'MW', 'Spin', Spin);

disp('Perform scalar directional harmonic-to-wavelet (lm2wav) transform with default parameters')
[f_wav, f_scal] = s2let_transform_analysis_lm2wav(flm, 'N', N, 'Upsample', true);
flm_rec = s2let_transform_synthesis_lm2wav(f_wav, f_scal, 'N', N, 'Upsample', true);
default = max(abs(flm-flm_rec))


%disp('Perform spin directional harmonic-to-wavelet (lm2wav) transform with default parameters')
%[f_wav, f_scal] = s2let_transform_analysis_lm2wav(flm, 'N', N, 'Spin', Spin, 'Upsample', true);
%flm_rec = s2let_transform_synthesis_lm2wav(f_wav, f_scal, 'N', N, 'Spin', Spin, 'Upsample', true);
%default = max(abs(flm-flm_rec))

disp('Perform spin directional harmonic-to-wavelet (lm2wav) transform with custom parameters')
[f_wav, f_scal] = s2let_transform_analysis_lm2wav(flm,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', false);
flm_rec = s2let_transform_synthesis_lm2wav(f_wav, f_scal,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', false);
default = max(abs(flm-flm_rec))

%disp('Perform scalar directional wavelet transform (analysis_mw) with default parameters')
%[f_wav, f_scal] = s2let_transform_analysis_mw(f, 'N', N, 'Upsample', true);
%f_rec = s2let_transform_synthesis_mw(f_wav, f_scal, 'N', N, 'Upsample', true);
%default = max(max(abs(f-f_rec)))

%disp('Perform spin directional wavelet transform with default parameters')
%[f_wav, f_scal] = s2let_transform_analysis_mw(f_s, 'N', N, 'Spin', Spin, 'Upsample', true);
%f_rec = s2let_transform_synthesis_mw(f_wav, f_scal, 'N', N, 'Spin', Spin, 'Upsample', true);
%default = max(max(abs(f_s-f_rec)))

disp('Perform spin directional wavelet transform (analysis and synthesis) with custom parameters')
[f_wav, f_scal] = s2let_transform_analysis_mw(f_s,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', false);
f_rec = s2let_transform_synthesis_mw(f_wav, f_scal,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', false);
default = max(max(abs(f_s-f_rec)))

%disp('Perform multiresolution scalar directional wavelet transform with default parameters')
%[f_wav, f_scal] = s2let_transform_analysis_mw(f, 'N', N);
%f_rec = s2let_transform_synthesis_mw(f_wav, f_scal, 'N', N);
%default = max(max(abs(f-f_rec)))

%disp('Perform multiresolution spin directional wavelet transform with default parameters')
%[f_wav, f_scal] = s2let_transform_analysis_mw(f_s, 'N', N, 'Spin', Spin);
%f_rec = s2let_transform_synthesis_mw(f_wav, f_scal, 'N', N, 'Spin', Spin);
%default = max(max(abs(f_s-f_rec)))

disp('Perform multiresolution spin directional wavelet transform with custom parameters')
[f_wav, f_scal] = s2let_transform_analysis_mw(f_s,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', true);
f_rec = s2let_transform_synthesis_mw(f_wav, f_scal,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', true);
default = max(max(abs(f_s-f_rec)))


disp('Constraint on flms to generate real signal')
for el = 0:L-1
   ind = el*el + el + 1;
   flm(ind,1) = real(flm(ind,1));
   for m = 1:el
      ind_pm = el*el + el + m + 1;
      ind_nm = el*el + el - m + 1;
      flm(ind_nm,1) = (-1)^m * conj(flm(ind_pm,1));
   end
end
disp('Construct the corresponding real signal on the sphere')
f_real = ssht_inverse(flm, L, 'Method', 'MW', 'Reality', true);


disp('Perform multiresolution real directional wavelet transform with default parameters')
[f_wav, f_scal] = s2let_transform_analysis_mw(f_real, 'N', N, 'Reality', true);
f_rec = s2let_transform_synthesis_mw(f_wav, f_scal, 'N', N, 'Reality', true);
default = max(max(abs(f_real-f_rec)))


%disp('Perform full resolution real directional wavelet transform with default parameters')
%[f_wav, f_scal] = s2let_transform_analysis_mw(f_real, 'N', N, 'Reality', true, 'Upsample', true);
%f_rec = s2let_transform_synthesis_mw(f_wav, f_scal, 'N', N, 'Reality', true, 'Upsample', true);
%default = max(max(abs(f_real-f_rec)))


% stop
% ---------------------
% Axisymmetric wavelets:
% ----------------------
disp('Generates band-limited function')
flm = zeros(L^2,1);
flm = rand(size(flm)) + sqrt(-1)*rand(size(flm));
flm = 2.*(flm - (1+sqrt(-1))./2);
disp('Construct the corresponding signal on the sphere')
f = ssht_inverse(flm, L, 'Method', 'MW');
disp('Construct the corresponding spin signal on the sphere')
f_s = ssht_inverse(flm, L, 'Method', 'MW', 'Spin', Spin);

disp('Perform axisym wavelet transform with default parameters')
[f_wav, f_scal] = s2let_transform_axisym_analysis_mw(f);
f_rec = s2let_transform_axisym_synthesis_mw(f_wav, f_scal);
default = max(max(abs(f-f_rec)))

disp('Perform axisym wavelet transform with multiresolution algorithm')
[f_wav, f_scal] = s2let_transform_axisym_analysis_mw(f, 'Upsample', false);
f_rec = s2let_transform_axisym_synthesis_mw(f_wav, f_scal, 'Upsample', false);
default_multires = max(max(abs(f-f_rec)))

disp('Perform axisym wavelet transform at full resolution')
[f_wav, f_scal] = s2let_transform_axisym_analysis_mw(f, 'Upsample', true);
f_rec = s2let_transform_axisym_synthesis_mw(f_wav, f_scal, 'Upsample', true);
default_fullres = max(max(abs(f-f_rec)))

disp('Perform axisym wavelet transform with custom parameters')
[f_wav, f_scal] = s2let_transform_axisym_analysis_mw(f, 'B', B, 'L', L, 'J_min', J_min);
f_rec = s2let_transform_axisym_synthesis_mw(f_wav, f_scal, 'B', B, 'L', L, 'J_min', J_min);
custom = max(max(abs(f-f_rec)))

disp('Constraint on flms to generate real signal')
for el = 0:L-1
   ind = el*el + el + 1;
   flm(ind,1) = real(flm(ind,1));
   for m = 1:el
      ind_pm = el*el + el + m + 1;
      ind_nm = el*el + el - m + 1;
      flm(ind_nm,1) = (-1)^m * conj(flm(ind_pm,1));
   end
end
disp('Construct the corresponding real signal on the sphere')
f_real = ssht_inverse(flm, L, 'Method', 'MW', 'Reality', true);


disp('Perform multiresolution real directional wavelet transform with default parameters')
[f_wav, f_scal] = s2let_transform_analysis_mw(f_real, 'N', N, 'Reality', true);
f_rec = s2let_transform_synthesis_mw(f_wav, f_scal, 'N', N, 'Reality', true);
default = max(max(abs(f_real-f_rec)))


disp('Perform real axisym wavelet transform with default parameters')
[f_wav_real, f_scal_real] = s2let_transform_axisym_analysis_mw(f_real, 'Reality', true);
f_real_rec = s2let_transform_axisym_synthesis_mw(f_wav_real, f_scal_real, 'Reality', true);
default = max(max(abs(f_real-f_real_rec)))

disp('Perform real axisym wavelet transform with multiresolution algorithm')
[f_wav_real, f_scal_real] = s2let_transform_axisym_analysis_mw(f_real, 'Upsample', false, 'Reality', true);
f_real_rec = s2let_transform_axisym_synthesis_mw(f_wav_real, f_scal_real, 'Upsample', false, 'Reality', true);
default_multires = max(max(abs(f_real-f_real_rec)))

disp('Perform real axisym wavelet transform at full resolution')
[f_wav_real, f_scal_real] = s2let_transform_axisym_analysis_mw(f_real, 'Upsample', true, 'Reality', true);
f_real_rec = s2let_transform_axisym_synthesis_mw(f_wav_real, f_scal_real, 'Upsample', true, 'Reality', true);
default_fullres = max(max(abs(f_real-f_real_rec)))

disp('Perform real axisym wavelet transform with custom parameters')
[f_wav_real, f_scal_real] = s2let_transform_axisym_analysis_mw(f_real, 'B', B, 'L', L, 'J_min', J_min, 'Reality', true);
f_real_rec = s2let_transform_axisym_synthesis_mw(f_wav_real, f_scal_real, 'B', B, 'L', L, 'J_min', J_min, 'Reality', true);
custom = max(max(abs(f_real-f_real_rec)))



