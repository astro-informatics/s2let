% s2let_curvelet_transform_test_lm2cur_cur2lm
%
% Run curvelet analysis (harmonic to curvelet space) and 
% synthesis (curvelet to harmonic space)
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

clear all;
close all;

% Curvelet parameters
Spin = 2;             % Spin value of curvelet 
L = 16;               % Angular band-limit
B = 2;                % B = 2 for dyadic sampling
J_min = 1;            % Minimum scale probed by curvelets
J =s2let_jmax(L, B);  % Maximum scale probed by curvelets =ceil(log L/ log B);

disp('Generates random band-limited function')
flm_gen = zeros(L^2,1);
flm_gen = rand(size(flm_gen)) + sqrt(-1)*rand(size(flm_gen));
flm_gen = 2.*(flm_gen - (1+sqrt(-1))./2);
disp('Construct the corresponding signal on the sphere')
f_gen = ssht_inverse(flm_gen, L,'Method', 'MW');
% Spin signals: 
flm_spin_gen= ssht_forward(f_gen, L, 'Spin', Spin, 'Method', 'MW');
f_spin_gen = ssht_inverse(flm_spin_gen, L, 'Spin', Spin, 'Method', 'MW');
disp('----------- ');


disp(' ');
disp('Curvelet transform: full resolution (Upsample: true): ');
% -----------------
% Signal analysis: (harmonic to curvelet space)
% -----------------
% N.B. 's2let_transform_curvelet_analysis_lm2cur.m' called 's2let_curvelet_tiling.m'
disp('Spin signal, Full resolution: analysis_lm2cur...')
[f_cur, f_scal] = s2let_transform_curvelet_analysis_lm2cur(flm_spin_gen,  ...
                                                           'B', B, 'L', L, ...
                                                           'J_min', J_min, ...
                                                           'Spin', Spin, ...
                                                           'Upsample', true); 

% -----------------
% Signal synthesis: (curvelet to harmonic space) 
% -----------------
disp('Spin signal, Full resolution: synthesis_cur2lm...')
flm_spin_rec= s2let_transform_curvelet_synthesis_cur2lm(f_cur, f_scal, ...
                                                        'B', B, 'L', L,...
                                                        'J_min', J_min, ...
                                                        'Spin', Spin, ...
                                                        'Upsample', true);


disp('Compute the re-constructed function via ssht_inverse ');
f_spin_rec = ssht_inverse(flm_spin_rec, L,'Spin', Spin,'Method', 'MW');
disp('- Test exact transform:');
disp('Check the difference between flm_gen and flm_rec:');
maxerr = max(abs(flm_spin_gen - flm_spin_rec))
disp('Check the difference between f_spin_gen and f_spin_rec: ');
maxerr = max(abs(f_spin_gen(:) - f_spin_rec(:)))
disp('----------- ');


% ================== MULTI-RESOLUTION ===================% 
disp(' Curvelet transform: multi-resolution (Upsample: false): ');
% -----------------
% Signal analysis: (harmonic to curvelet space)
% -----------------

[f_cur, f_scal] = s2let_transform_curvelet_analysis_lm2cur(flm_spin_gen,  ...
                                                           'B', B, 'L', L, ...
                                                           'J_min', J_min, ...
                                                           'Spin', Spin, ...
                                                           'Sampling', 'MW',...
                                                           'Upsample', false);

% -----------------
% Signal synthesis: (curvelet to harmonic space)
% -----------------
flm_spin_rec= s2let_transform_curvelet_synthesis_cur2lm(f_cur, f_scal, ...
                                                        'B', B, 'L', L,...
                                                        'J_min', J_min, ...
                                                        'Spin', Spin, ...
                                                        'Sampling', 'MW', ...
                                                        'Upsample', false);


disp('Compute the re-constructed function via ssht_inverse ');
f_spin_rec = ssht_inverse(flm_spin_rec, L, 'Spin', Spin, 'Method', 'MW');

disp('- Test exact transform:');
disp('Check the difference between flm_gen and flm_rec:');
maxerr = max(abs(flm_spin_gen - flm_spin_rec))
disp('Check the difference between f_spin_gen and f_spin_rec: ');
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
reality = true; 


% ================== FULL-RESOLUTION ===================% 
disp(' Curvelet transform: Full-resolution (Upsample: true): ');
% -----------------
% Signal analysis: (harmonic to curvelet space) 
% -----------------
disp('Real signal, Full resolution: analysis_lm2cur...')
upsample = true; 
[f_cur, f_scal]= s2let_transform_curvelet_analysis_lm2cur(flm_gen, ...
                                                          'B', B, 'L', L,...
                                                          'J_min', J_min,...
                                                          'Reality', reality, ...
                                                          'Upsample', upsample);
% -----------------
% Signal synthesis: (curvelet to harmonic space)
% -----------------
disp('Real signal, Full resolution: synthesis_cur2lm...')
flm_rec  = s2let_transform_curvelet_synthesis_cur2lm(f_cur, f_scal,  ...
                                                     'B', B, 'L', L, ...
                                                     'J_min', J_min,...
                                                     'Reality', reality,...
                                                     'Upsample', upsample);
f_gen = ssht_inverse(flm_gen, L,'Method', 'MW');
f_rec = ssht_inverse(flm_rec, L,'Method', 'MW');
disp('- Test exact transform:');
disp('Check the difference between flm_gen and flm_rec:');
maxerr = max(abs(flm_gen - flm_rec))
disp('Check the difference between f_gen and f_rec: ');
maxerr = max(abs(f_gen(:) - f_rec(:)))
disp('----------- ');


% ================== MULTI-RESOLUTION ===================% 
disp(' Curvelet transform: Multi-resolution (Upsample: false): ');
% -----------------
% Signal analysis: (harmonic to curvelet space)
% -----------------
upsample = false; 
disp('Real signal, Multi-reolution: analysis_lm2cur...')
[f_cur, f_scal]= s2let_transform_curvelet_analysis_lm2cur(flm_gen, ...
                                                          'B', B, 'L', L, ...
                                                          'J_min', J_min,...
                                                          'Reality', reality, ...
                                                          'Upsample', upsample);
% -----------------
% Signal synthesis: (curvelet to harmonic space)
% -----------------
disp('Real signal, Multi-reolution: synthesis_cur2lm...')
flm_rec  = s2let_transform_curvelet_synthesis_cur2lm(f_cur, f_scal,  ...
                                                     'B', B, 'L', L, ...
                                                     'J_min', J_min,...
                                                     'Reality', reality,...
                                                     'Upsample',upsample);
f_gen = ssht_inverse(flm_gen, L,'Method', 'MW');
f_rec = ssht_inverse(flm_rec, L,'Method', 'MW');
disp('- Test exact transform:');
disp('Check the difference between flm_gen and flm_rec:');
maxerr = max(abs(flm_gen - flm_rec))
disp('Check the difference between f_gen and f_rec: ');
maxerr = max(abs(f_gen(:) - f_rec(:)))
disp('----------- ');