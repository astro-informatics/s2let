

% ---------------
% Define curvelet parameters: 
% ---------------
B = 2;
L = 16; %32; 
Spin = 0; 
N= L;     % Since m=l, the azimuthal band limit N = overall band limit L
J_min = 3; % minimum and maximum scale probed by wavelets 
J =s2let_jmax(L, B);  %=ceil(log L/ log B) =5; 

f_scal = zeros(L, 2*L-1);

f_cur=cell(J+1-J_min);
for j = J_min:J  
band_limit = min([ s2let_bandlimit(j,J_min,B,L) L ]);  
Nj = band_limit
f_cur{j-J_min+1} = zeros(2*Nj-1, L, 2*L-1);
end 
% for j= j_min
j = randi(J-J_min+1)+J_min-1;
band_limit = min([ s2let_bandlimit(j,J_min,B,L) L ]);  
Nj = band_limit

gamma= randi(2*Nj)-1; 
alpha = 1; %randi(L+1)-1; 
beta = 1;%  randi(2*L)-1; 
f_cur{j-J_min+1}(gamma,alpha,beta) = 1.;

f_rec = s2let_transform_spin_curvelet_synthesis_cur2px(f_cur, f_scal, ...
                                                      'B', B, 'L', L, ...
                                                      'J_min', J_min, ...
                                                      'Spin', Spin, ...
                                                      'Reality', false, ...
                                                      'Upsample', true, ...
                                                      'SpinLowered', false, ...
                                                      'SpinLoweredFrom', 0,...
                                                 	  'Sampling', 'MW');

                                                  
figure
ssht_plot_mollweide(real(f_rec),L, 'Mode', 1);
%ssht_plot_sphere(real(f_rec), L, 'Type', 'colour', 'Lighting', true);
                                                  