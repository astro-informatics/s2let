function  [cur_lm scal_l] = s2let_curvelet_tiling(B, L, J_min, varargin)

% s2let_wavelet_tiling - Compute tiling in harmonic space.
% -- CURVELETS on the sphere.
%
% Default usage :
%
%   [cur_lm scal_l] = s2let_spin0curvelet_tiling(B, L, J_min, <options>)
%
% cur_lm is an array containing the curvelets spherical harmonic coefficients.
% scal_l is an array containing the scaling function spherical harmonic coefficients (l only).
% B is the wavelet parameter,
% L is the harmonic band-limit;
% J_min the first wavelet to be used.
%  
% % Valid options include:
%
%  'Spin'        = { Spin; (default=0) }
%  'SpinLowered' = { true  [Apply normalisation factors for spin-lowered
%                           wavelets and scaling function.],
%                    false [Apply the usual normalisation factors such
%                           that the wavelets fulfil the admissibility
%                           condition (default)]}
%  'SpinLoweredFrom' = [integer; if the SpinLowered option is used, this
%                       option indicates which spin number the wavelets
%                       should be lowered from (default = 0)]
%
% -----------------------------------------------------------
% Log: 
% -  constructed by Jennifer Y H Chan on 5th June 2015  
% -----------------------------------------------------------
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

p = inputParser;
p.addRequired('B', @isnumeric);
p.addRequired('L',  @isnumeric);
p.addRequired('J_min', @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('SpinLowered', false, @islogical);
p.addParamValue('SpinLoweredFrom', 0, @isnumeric);
p.parse(B, L, J_min, varargin{:});
args = p.Results;

% N = L; 
J = s2let_jmax(L, B);
Spin = args.Spin;
original_spin= args.SpinLoweredFrom; 

% Effectively ignore original_spin if we don't use spin-lowered wavelets.
if (args.SpinLowered~= 0)
original_spin = 0 ;  
end 
el_min = max(abs(Spin), abs(original_spin));

% ----------
% Curvelet directional component 
% ----------
signs = zeros(L,1);  %, size('double')
s_lm = zeros(L^2,1); 
% Perform precomputation.
for (m=1:2:L-1) 
  signs(m)   = -1.0;
  signs(m+1) = 1.0;
end 
% Skip the s_00 component as it is zero 
for el = 1:L-1 
% for el = Spin*Spin+1:L-1 
m = el;  % automatically satisfy the condition : m < N 
% for positive m
ind_pm = ssht_elm2ind(el, m);
% Curvelet coefficients:
s_lm(ind_pm)= sqrt(1./2.);
% for negative m
ind_nm = ssht_elm2ind(el, -m);
%s_lm(ind_nm)= (-1)^m*conj(s_lm(ind_pm));
s_lm(ind_nm)= signs(m)*conj(s_lm(ind_pm));
end


% `````````````
% Check the zero mean condition for the curvelet directional components 
% `````````````
%{
ind = 1; 
error_on_slm_tiling = 0; 
% % Skip the s_00 component as it is zero 
for el = 0:L-1
sum = 0.; 
for m= -el:el
 sum = sum + s_lm(ind)*conj(s_lm(ind));
 ind = ind +1; 
end
end
error_on_slm_tiling = error_on_slm_tiling +sum - 1.0
%} 


% ----------
% Curvelet angular components:
% ----------
[kappa kappa0] =  s2let_transform_axisym_tiling(B, L, J_min);
ind_pm = 0;
ind_nm = 0;
%ind_pm = el_min*el_min;
%ind_pm = el_min*el_min;
for j = J_min:J
% for el = 1:L-1  %s_00 is zero
for el = el_min+1:L-1 
 m = el;
 % for positive m
 ind_pm = ssht_elm2ind(el, m);
 % Curvelet coefficients:
 % s_lm= sqrt(1./2.);
 % cur_lm{j-J_min+1}(ind_pm) =sqrt(1./2.) * sqrt((2*el+1)/(8.0*pi*pi))* kappa(j+1,el+1);
 cur_lm{j-J_min+1}(ind_pm) = s_lm(ind_pm) * sqrt((2*el+1)/(8.0*pi*pi))* kappa(j+1,el+1);
 % for negative m
 ind_nm = ssht_elm2ind(el, -m);
 cur_lm{j-J_min+1}(ind_nm) =  signs(m) * conj(cur_lm{j-J_min+1}(ind_pm));
 % cur_lm{j-J_min+1}(ind_nm) =(-1)^m* conj(cur_lm{j-J_min+1}(ind_pm));
  if (args.SpinLowered~= 0)  
    s2let_spin_lowered_norm_factor = s2let_spin_lowered_normalization(el, 'original_spin',original_spin);
    cur_lm{j-J_min+1}(ind_pm) = cur_lm{j-J_min+1}(ind_pm)+ s2let_spin_lowered_norm_factor ;
    cur_lm{j-J_min+1}(ind_nm) = cur_lm{j-J_min+1}(ind_nm)+ s2let_spin_lowered_norm_factor ;
  end 
 end
end 

% ----------
% Scaling Function
% ----------
scal_l = zeros(L^2,1);
for el = 0:L-1
%for l = el_min:L-1
 scal_l(el^2+el+1,1) = sqrt((2*el+1)/(4.0*pi)) *kappa0(el+1);
 if (args.SpinLowered~= 0)  
  s2let_spin_lowered_norm_factor = s2let_spin_lowered_normalization(el, 'original_spin',original_spin);
  scal_l(el^2+el+1,1) = scal_l(el^2+el+1,1) *s2let_spin_lowered_norm_factor ;
 end 
end


% `````````````
% Check admissibility condition (OR call external function
% 's2let_check_cur_tiling.m', run 's2let_fulltest.m')
% `````````````
%{
Spin = 0; 
% Scaling function
identity = zeros(1,L);
for el=abs(Spin):L-1
	% identity(1,el+1) = identity(1,el+1) + 4*pi/(2*el+1) * kappa0_cur(el+1,1) * conj(kappa0_cur(el+1,1));
    identity(1,el+1) = identity(1,el+1)+(4.*pi/(2*el+1))*scal_l(el^2+el+1,1)*conj(scal_l(el^2+el+1,1)); 
end
% Curvelet functions
for j= 0:J %J_min: J 
	ind = Spin*Spin + 1;
	for el=abs(Spin):L-1
		for m= -el:el
           identity(1,el+1) = identity(1,el+1)+(8.*pi^2/(2*el+1))* cur_lm{j-J_min+1}(1,ind)*conj(cur_lm{j-J_min+1}(1,ind));
		   ind = ind + 1;
		end
	end
end
error_on_cur_tiling = 0;
for el=abs(Spin):L-1
    error_on_cur_tiling = error_on_cur_tiling + identity(1,el+1) - 1.0;
end
%}

end
