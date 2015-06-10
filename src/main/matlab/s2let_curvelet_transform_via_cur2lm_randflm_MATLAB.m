% s2let_curvelet_matlab_prototype
% Run curvelet analysis and synthesis 
% of randomly generated signals f 

clear all;
close all;

% Main parameters
Spin = 0;
L = 16;
N = L; % For curvelets: N=L  %N=1 for axisymmetric curvelets
B = 2;
J_min = 1;  %N.B. J_min = 2 : NOT EXACT
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


% -----------------
% Signal analysis: (harmonic to pixel space) 
% -----------------
% Call matlab function analysis_lm2cur
[f_cur, f_scal] = s2let_transform_curvelet_analysis_lm2cur(flm_gen,  ...
                                                  'B', B, 'L', L, ...
                                                  'J_min', J_min, 'Spin', Spin, ...
                                                  'Sampling', 'MW',...
                                                  'Upsample', true);
                                              
% -----------------
% Signal synthesis: (pixel to harmonic space) 
% -----------------
% Call s2let_transform_synthesis_lmn2lm
flm_rec= s2let_transform_curvelet_synthesis_cur2lm(f_cur, f_scal, ...
                                           'B', B, 'L', L,...
                                           'J_min', J_min, ...
                                           'Spin', Spin, ...
                                           'Sampling', 'MW', ...
                                           'Upsample', true);


                                       
disp('Compute the re-constructed function via ssht_inverse ');
f_rec = ssht_inverse(flm_rec, L, 'Method', 'MW'); 

disp('Both analysis and synthesis are done! ');
disp('');
disp('- Test exact transform: check the difference between flm_gen and flm_rec:');
maxerr = max(abs(flm_gen - flm_rec))
disp('Check the difference between f_gen and f_rec: ');
maxerr = max(abs(f_gen(:) - f_rec(:)))
