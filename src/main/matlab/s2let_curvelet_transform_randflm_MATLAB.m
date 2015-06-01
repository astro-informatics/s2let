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

disp('Generates random band-limited function')
flm_gen = zeros(L^2,1);
flm_gen = rand(size(flm_gen)) + sqrt(-1)*rand(size(flm_gen));
flm_gen = 2.*(flm_gen - (1+sqrt(-1))./2);

disp('Construct the corresponding signal on the sphere')
f_gen = ssht_inverse(flm_gen, L, 'Method', 'MW');


% ---------------
% Tile curvelets:
% ---------------
%***** step 1a) call curvelet- and scaling-function- generating functions 
disp('Tile curvelets in harmonic space (cur_lm, scal_l)')
[psi_lm phi_l] = s2let_curvelet_tiling(B, L, N, Spin, J_min);
for j = J_min:J, 
 cur_lm{j-J_min+1} = psi_lm(:,j+1);  
end
%***** step 1b) compute the scaling coefficients (no j-dependence except on J_min)
scal_l = zeros(L^2,1);
for l = 0:L-1,
scal_l(l^2+l+1,1) = phi_l(l+1);
end


% -----------------
% Signal analysis: (harmonic to pixel space) 
% -----------------
% Call matlab function analysis_lm2cur
[f_cur, f_scal] = s2let_transform_analysis_lm2cur(flm_gen,  ...
                                                  'B', B, 'L', L, ...
                                                  'Spin', Spin, 'J_min', J_min, ...
                                                  'N', N, 'Upsample', false, 'Sampling', 'MW');
                                              
% -----------------
% Signal synthesis: (pixel to harmonic space) 
% -----------------
% Call s2let_transform_synthesis_cur2lm
% Call s2let_transform_synthesis_lmn2lm
[flm_cur_syn, flm_scal_syn] = s2let_transform_synthesis_lmn2lm(f_cur_lmn, f_scal_lm, ...
                                                               cur_lm, scal_l, flm_gen, ...
                                                               'B', B, 'L', L, 'J_min', J_min, ...
                                                               'N', N, 'Upsample', false, 'Spin', 0);

disp('Compute the re-constructed function via ssht_inverse ');
f_rec = ssht_inverse(flm_rec, L, 'Method', 'MW'); 

disp('Both analysis and synthesis are done! ');
disp('');
disp('- Test exact transform: check the difference between flm_gen and flm_rec:');
maxerr = max(abs(flm_gen - flm_rec))
%disp('Check the difference between f_gen and f_rec: ');
%maxerr = max(abs(f_gen - f_rec))
