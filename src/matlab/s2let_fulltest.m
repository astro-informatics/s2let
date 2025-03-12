% s2let_fulltest
% Run all exactness tests for the MW sampling,
% all wavelet transforms must reconstruct the input maps
% at floating-point precision. Various parameters are tested.
%
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

clear all;
close all;

% Main parameters
L = 32;
B = 4;
N = 4;
Spin = 0;
J_min = 1;
J = s2let_jmax(L, B);

disp('Checks tiling of harmonic space for axysimmetric wavelets')
[kappa kappa0] = s2let_transform_axisym_tiling(B, L, J_min);
error_on_axisym_tiling = s2let_check_axisym_tiling(kappa, kappa0, L, J)

disp('Checks tiling of harmonic space for directional wavelets')
[psi phi] = s2let_wavelet_tiling(B, L, N, Spin, J_min);
error_on_axisym_tiling = s2let_check_tiling(psi, phi, L, Spin, J)

disp('Generates band-limited function')
flm = zeros(L^2,1);
flm = rand(size(flm)) + sqrt(-1)*rand(size(flm));
flm = 2.*(flm - (1+sqrt(-1))./2);
disp('Construct the corresponding signal on the sphere')
f = ssht_inverse(flm, L, 'Method', 'MW');
disp('Construct the corresponding spin signal on the sphere')
f_s = ssht_inverse(flm, L, 'Method', 'MW', 'Spin', Spin);

% disp('Perform scalar directional harmonic-to-wavelet transform with default parameters')
% [f_wav, f_scal] = s2let_transform_analysis_lm2wav(flm, 'N', N, 'Upsample', true);
% flm_rec = s2let_transform_synthesis_lm2wav(f_wav, f_scal, 'N', N, 'Upsample', true);
% default = max(abs(flm-flm_rec))

% disp('Perform spin directional harmonic-to-wavelet transform with default parameters')
% [f_wav, f_scal] = s2let_transform_analysis_lm2wav(flm, 'N', N, 'Spin', Spin, 'Upsample', true);
% flm_rec = s2let_transform_synthesis_lm2wav(f_wav, f_scal, 'N', N, 'Spin', Spin, 'Upsample', true);
% default = max(abs(flm-flm_rec))

% disp('Perform spin directional harmonic-to-wavelet transform with custom parameters')
% [f_wav, f_scal] = s2let_transform_analysis_lm2wav(flm,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', false);
% flm_rec = s2let_transform_synthesis_lm2wav(f_wav, f_scal,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', false);
% default = max(abs(flm-flm_rec))

disp('Perform scalar directional wavelet transform with default parameters')
[f_wav, f_scal] = s2let_transform_analysis_mw(f, 'N', N, 'Upsample', true);
f_rec = s2let_transform_synthesis_mw(f_wav, f_scal, 'N', N, 'Upsample', true);
default = max(max(abs(f-f_rec)))

disp('Perform spin directional wavelet transform with default parameters')
[f_wav, f_scal] = s2let_transform_analysis_mw(f_s, 'N', N, 'Spin', Spin, 'Upsample', true);
f_rec = s2let_transform_synthesis_mw(f_wav, f_scal, 'N', N, 'Spin', Spin, 'Upsample', true);
default = max(max(abs(f_s-f_rec)))

disp('Perform spin directional wavelet transform with custom parameters')
[f_wav, f_scal] = s2let_transform_analysis_mw(f_s,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', false);
f_rec = s2let_transform_synthesis_mw(f_wav, f_scal,  'B', B, 'L', L, 'J_min', J_min, 'N', N, 'Spin', Spin, 'Upsample', false);
default = max(max(abs(f_s-f_rec)))

disp('Perform multiresolution scalar directional wavelet transform with default parameters')
[f_wav, f_scal] = s2let_transform_analysis_mw(f, 'N', N);
f_rec = s2let_transform_synthesis_mw(f_wav, f_scal, 'N', N);
default = max(max(abs(f-f_rec)))

disp('Perform multiresolution spin directional wavelet transform with default parameters')
[f_wav, f_scal] = s2let_transform_analysis_mw(f_s, 'N', N, 'Spin', Spin);
f_rec = s2let_transform_synthesis_mw(f_wav, f_scal, 'N', N, 'Spin', Spin);
default = max(max(abs(f_s-f_rec)))

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


disp('Perform full resolution real directional wavelet transform with default parameters')
[f_wav, f_scal] = s2let_transform_analysis_mw(f_real, 'N', N, 'Reality', true, 'Upsample', true);
f_rec = s2let_transform_synthesis_mw(f_wav, f_scal, 'N', N, 'Reality', true, 'Upsample', true);
default = max(max(abs(f_real-f_rec)))

stop


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

