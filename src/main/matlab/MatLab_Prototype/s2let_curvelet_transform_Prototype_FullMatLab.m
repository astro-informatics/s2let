% s2let_curvelet_transform_via_lmn2lm_randflm_MATLAB
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
flm_init = zeros(L^2,1);
flm_init = rand(size(flm_init)) + sqrt(-1)*rand(size(flm_init));
flm_init = 2.*(flm_init - (1+sqrt(-1))./2);

disp('Construct the corresponding signal on the sphere')
f_init = ssht_inverse(flm_init, L, 'Method', 'MW');
%}

flm_init = zeros(L^2,1);
%
disp('read from file flm_init'); 
fid= fopen('/Users/jenniferyhchan/WaveletsCode_PhD/s2let_curvelets_MATLAB/1_cur_flm_randgen_mw_test.dat');
rawData=fscanf(fid, '%f, %f',[2 256]);
fclose(fid);
complexData=complex(rawData(1,:),rawData(2,:));
flm_init= complexData.' ;  % Non-conjugate transpose
%
disp('Construct the corresponding signal via ssht_inverse'); 
f_init = ssht_inverse(flm_init, L, 'Method', 'MW');
%
disp('Construct the corresponding flm of the signal via ssht_forward'); 
flm_init= ssht_forward(f_init, L, 'Method', 'MW');

% ---------------
% Tile curvelets:
% ---------------
% Call curvelet- and scaling-function- generating functions 
disp('curvelet_tiling: Tile curvelets in harmonic space (cur_lm, scal_l)')
[cur_lm scal_l] = s2let_curvelet_tiling(B, L, J_min, ...
                                        'Spin', Spin, 'SpinLowered', false,  'SpinLoweredFrom', 0);

% -----------------
% Signal analysis: (harmonic to Wigner space)
% -----------------
%{
% Call matlab function analysis_lm2lmn
[f_cur_lmn, f_scal_lm] = s2let_transform_analysis_lm2lmn(flm_init, cur_lm, scal_l, ...
                                                         'B', B, 'L', L, 'J_min', J_min,...
                                                         'Spin', Spin,  'Reality', false, 'Upsample', false,...
                                                         'SpinLowered', false,  'SpinLoweredFrom', 0, ...
                                                         'Sampling', 'MW');
%}
%
% Curvelet Contribution
%
for j = J_min:J,
 f_cur_lmn{j-J_min+1} = zeros((2*N-1)*L*L,1);
 ind_ln=0;
 ind = 0;
 ind_lmn = 0;
 for n = -N+1:N-1,
    for el = abs(n):L-1,
        ind_ln = ssht_elm2ind(el, n);
        psi = 8.*pi*pi/(2.*el+1) *conj(cur_lm{j-J_min+1}(ind_ln));
        for m = -el:el,
        ind = ssht_elm2ind(el, m); 
        ind_lmn = so3_elmn2ind(el,m,n,L,N);
        f_cur_lmn{j-J_min+1}(ind_lmn) =  flm_init(ind) * psi;
        end
    end
 end
end
%
% Scaling function Contribution
%
disp('Compute scaling function f_scal_lm=flm_init(lm_ind) * phi ')
%f_scal_lm = zeros(L^2,1)
lm_ind=0;
for el = 0:L-1,
 phi = sqrt(4.0*pi/(2.*el+1)) * scal_l(el^2+el+1,1);  %
 for m = -el:el,
 lm_ind=ssht_elm2ind(el, m);
 f_scal_lm(lm_ind) = flm_init(lm_ind) * phi;
 end
end


%
disp('so3_inverse and so3_forward')
%
for j = J_min:J,
f_cur{j-J_min+1} = so3_inverse(f_cur_lmn{j-J_min+1}, L, N, 'Sampling', 'MW');
f_cur_lmn{j-J_min+1} = so3_forward(f_cur{j-J_min+1}, L, N,'Sampling', 'MW');
end  
f_scal= ssht_inverse(f_scal_lm, L,'Method', 'MW');
f_scal_lm = ssht_forward(f_scal, L,'Method', 'MW');
                                            
% -----------------
% Signal synthesis: (Wigner to harmonic space)
% -----------------
%{
% Call s2let_transform_synthesis_lmn2lm
[flm_cur_syn, flm_scal_syn] = s2let_transform_synthesis_lmn2lm(f_cur_lmn, f_scal_lm, ...
                                                               cur_lm, scal_l, flm_init, ...
                                                               'B', B, 'L', L, 'J_min', J_min,...
                                                               'Spin', Spin,  'Reality', false, 'Upsample', false,...
                                                               'SpinLowered', false,  'SpinLoweredFrom', 0, ...
                                                               'Sampling', 'MW');
%}
disp('Signal synthesis: ');
disp('Compute flm_cur_syn from f_cur_lmn');
flm_cur_syn=zeros(L^2,1);
for j = J_min:J,
 ind_ln =0;
 ind=0;
 ind_lmn=0;
 for en = -N+1:N-1,
    for el = abs(en):L-1,
        ind_ln = ssht_elm2ind(el, en);
        psi = (cur_lm{j-J_min+1}(ind_ln));
        for m = -el:el,
        ind = ssht_elm2ind(el, m);
        ind_lmn = so3_elmn2ind(el,m,en,L,N);
        flm_cur_syn(ind) =flm_cur_syn(ind)+ f_cur_lmn{j-J_min+1}(ind_lmn)* psi;
        end 
    end
 end
end
disp('Compute flm_scal_syn ')
%% For scaling function - sum over m: 
flm_scal_syn=zeros(L^2,1);
lm_ind=0;
 for el = 0:L-1,
  phi = sqrt(4.*pi/(2.*el+1))*scal_l(el^2+el+1,1);
  for m = -el:el, 
  lm_ind=ssht_elm2ind(el, m);
  flm_scal_syn(lm_ind) =  flm_scal_syn(lm_ind)+ f_scal_lm(lm_ind)* phi;
  end
 end
% f_scal_l=sum(flm_test_scal)

disp('Sum: flm_cur_syn+flm_scal_syn ');
flm_rec = flm_scal_syn + flm_cur_syn;

                                       
disp('Compute the re-constructed function via ssht_inverse ');
f_rec = ssht_inverse(flm_rec, L, 'Method', 'MW'); 

disp('Both analysis and synthesis are done! ');
disp('');
disp('- Test exact transform: check the difference between flm_init and flm_rec:');
maxerr = max(abs(flm_init - flm_rec))
disp('Check the difference between f_init and f_rec: ');
maxerr = max(abs(f_init(:) - f_rec(:)))
