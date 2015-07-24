% s2let_curvelet_matlab_prototype
% Run curvelet analysis and synthesis 
% of randomly generated signals f 

clear; % all;
% close all;

% Main parameters
Spin = 0;
L = 16;
N = L; % For curvelets: N=L  %N=1 for axisymmetric curvelets
B = 2;
J_min = 1;  %N.B. J_min = 2 : NOT EXACT
J =s2let_jmax(L, B);  %=ceil(log L/ log B);  

disp('Generates random band-limited function')
flm_gen = zeros(L^2,1);
flm_gen = rand(size(flm_gen)) + sqrt(-1)*rand(size(flm_gen));
flm_gen = 2.*(flm_gen - (1+sqrt(-1))./2);
disp('Construct the corresponding signal on the sphere')
f_gen = ssht_inverse(flm_gen, L,'Method', 'MW');
f_spin_gen = ssht_inverse(flm_gen, L, 'Spin', Spin, 'Method', 'MW');
flm_spin_gen= ssht_forward(f_gen, L, 'Spin', Spin, 'Method', 'MW');
%}

%{
flm_gen = zeros(L^2,1);
disp('read from file flm_gen'); 
fid= fopen('/Users/jenniferyhchan/WaveletsCode_PhD/s2let_curvelets_MATLAB/1_cur_flm_randgen_mw_test.dat');
rawData=fscanf(fid, '%f, %f',[2 256]);
fclose(fid);
complexData=complex(rawData(1,:),rawData(2,:));
flm_gen= complexData.' ;  % Non-conjugate transpose
disp('Construct the corresponding signal via ssht_inverse'); 
f_spin_gen = ssht_inverse(flm_gen, L,  'Spin', Spin, 'Method', 'MW');
flm_spin_gen= ssht_forward(f_spin_gen, L,  'Spin', Spin,'Method', 'MW');
%}

%{
disp(' ');
disp('----------- ');
disp('Curvelet transform: full resolution (Upsample: true): ');
% -----------------
% Signal analysis: (harmonic to pixel space)
% -----------------
% N.B. 's2let_transform_curvelet_analysis_lm2cur.m'
% called 's2let_curvelet_tiling.m'
[f_cur, f_scal] = s2let_transform_spin_curvelet_analysis_lm2cur(flm_spin_gen,  ...
                                                           'B', B, 'L', L, ...
                                                           'J_min', J_min, ...
                                                           'Spin', Spin, ...
                                                           'Sampling', 'MW',...
                                                           'Upsample', true);

% -----------------
% Signal synthesis: (pixel to harmonic space)
% -----------------
flm_spin_rec= s2let_transform_spin_curvelet_synthesis_cur2lm(f_cur, f_scal, ...
                                                   'B', B, 'L', L,...
                                                   'J_min', J_min, ...
                                                   'Spin', Spin, ...
                                                   'Sampling', 'MW', ...
                                                   'Upsample', true);


disp('Compute the re-constructed function via ssht_inverse ');
f_spin_rec = ssht_inverse(flm_spin_rec, L,'Spin', Spin,'Method', 'MW');

disp('- Test exact transform:');
disp('Check the difference between flm_gen and flm_rec:');
maxerr = max(abs(flm_spin_gen - flm_spin_rec))
disp('Check the difference between f_gen and f_rec: ');
maxerr = max(abs(f_spin_gen(:) - f_spin_rec(:)))



% ================== MULTI-RESOLUTION ===================% 
disp('----------- ');
disp(' Curvelet transform: multi-resolution (Upsample: false): ');
% -----------------
% Signal analysis: (harmonic to pixel space)
% -----------------
% N.B. 's2let_transform_spin_curvelet_analysis_lm2cur.m'
% called 's2let_spin_curvelet_tiling.m'
[f_cur, f_scal] = s2let_transform_spin_curvelet_analysis_lm2cur(flm_spin_gen,  ...
                                                                'B', B, 'L', L, ...
                                                                'J_min', J_min, ...
                                                                'Spin', Spin, ...
                                                                'Sampling', 'MW',...
                                                                'Upsample', false);

% -----------------
% Signal synthesis: (pixel to harmonic space)
% -----------------
flm_spin_rec= s2let_transform_spin_curvelet_synthesis_cur2lm(f_cur, f_scal, ...
                                                   'B', B, 'L', L,...
                                                   'J_min', J_min, ...
                                                   'Sampling', 'MW', ...
                                                   'Upsample', false);


disp('Compute the re-constructed function via ssht_inverse ');
f_spin_rec = ssht_inverse(flm_spin_rec, L, 'Spin', Spin, 'Method', 'MW');

disp('- Test exact transform:');
disp('Check the difference between flm_gen and flm_rec:');
maxerr = max(abs(flm_spin_gen - flm_spin_rec))
disp('Check the difference between f_gen and f_rec: ');
maxerr = max(abs(f_spin_gen(:) - f_spin_rec(:)))
disp('----------- ');

%}


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
% isreal(flm_gen)  %0 i.e. No

% -----------------
% Signal analysis: (harmonic to wigner space) - Full resolution 
% -----------------
reality = true; 
upsample = true; 
% Call matlab function 's2let_transform_spin_curvelet_analysis_lm2cur'
disp('analysis_lm2cur...')
[f_cur, f_scal]= s2let_transform_spin_curvelet_analysis_lm2cur(flm_gen, ...
                                                               'B', B, 'L', L, 'J_min', J_min,...
                                                               'Reality', reality, ...
                                                               'Upsample', upsample,...
                                                               'Sampling', 'MW');
% -----------------
% Signal synthesis: (pixel to harmonic space) - Full resolution 
% -----------------
% Call matlab function 's2let_transform_spin_curvelet_synthesis_cur2lm'
disp('synthesis_cur2lm...')
flm_rec  = s2let_transform_spin_curvelet_synthesis_cur2lm(f_cur, f_scal,  ...
                                                          'B', B, 'L', L, 'J_min', J_min,...
                                                          'Reality', reality,...
                                                          'Upsample', upsample,...
                                                          'Sampling', 'MW');
f_gen = ssht_inverse(flm_gen, L,'Method', 'MW');
f_rec = ssht_inverse(flm_rec, L,'Method', 'MW');
disp('- Test exact transform:');
disp('Check the difference between flm_gen and flm_rec:');
maxerr = max(abs(flm_gen - flm_rec))
disp('Check the difference between f_gen and f_rec: ');
maxerr = max(abs(f_gen(:) - f_rec(:)))
disp('----------- ');
% -----------------
% Signal analysis: (harmonic to wigner space) - Multi-resolution 
% -----------------
upsample = false; 
disp('Multiresolution...')
% Call matlab function 's2let_transform_spin_curvelet_analysis_lm2cur'
disp('analysis_lm2cur...')
[f_cur, f_scal]= s2let_transform_spin_curvelet_analysis_lm2cur(flm_gen, ...
                                                                'B', B, 'L', L, 'J_min', J_min,...
                                                                'Reality', reality, ...
                                                                'Upsample', upsample,...
                                                                'Sampling', 'MW');
% -----------------
% Signal synthesis: (pixel to harmonic space)
% -----------------
% Call matlab function 's2let_transform_spin_curvelet_synthesis_cur2lm'
disp('synthesis_cur2lm...')
flm_rec  = s2let_transform_spin_curvelet_synthesis_cur2lm(f_cur, f_scal,  ...
                                                          'B', B, 'L', L, 'J_min', J_min,...
                                                          'Reality', reality,...
                                                          'Upsample',upsample,...
                                                          'Sampling', 'MW');
f_gen = ssht_inverse(flm_gen, L,'Method', 'MW');
f_rec = ssht_inverse(flm_rec, L,'Method', 'MW');
disp('- Test exact transform:');
disp('Check the difference between flm_gen and flm_rec:');
maxerr = max(abs(flm_gen - flm_rec))
disp('Check the difference between f_gen and f_rec: ');
maxerr = max(abs(f_gen(:) - f_rec(:)))
disp('----------- ');