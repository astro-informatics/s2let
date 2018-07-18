function  [cur_lm scal_l] = s2let_spin_curvelet_tiling(B, L, J_min, varargin)
% s2let_spin_curvelet_tiling - Tile scaling functions and curvelets in harmonic space.
%
% Default usage :
%
%   [cur_lm scal_l] = s2let_spin_curvelet_tiling(B, L, J_min, <options>)
%
% cur_lm is an array containing the unrotated curvelets spherical harmonic coefficients.
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
% S2LET package to perform wavelets transform on the Sphere.
% Copyright (C) 2012-2016  Boris Leistedt, Jennifer Chan & Jason McEwen
% See LICENSE.txt for license details
% -----------------------------------------------------------

p = inputParser;
p.addRequired('B', @isnumeric);
p.addRequired('L',  @isnumeric);
p.addRequired('J_min', @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('SpinLowered', false, @islogical);
p.addParamValue('SpinLoweredFrom', 0, @isnumeric);
p.parse(B, L, J_min, varargin{:});
args = p.Results;

J = s2let_jmax(L, B);
Spin = args.Spin; 
if (args.SpinLowered ~= 0) 
    original_spin= args.SpinLoweredFrom;  % For spin-lowered curvelet: (i.e. use scalar curvelets for the transform : spin =0, SpinLoweredFrom = e.g. 2)
else 
    original_spin = 0 ;  % (default - not to use spin-lowered wavelets).    
end 

% ----------
% Curvelet directional component s_lm 
% ----------
signs = zeros(L,1); 
s_lm = zeros(L^2,1); 
% Perform precomputation.
for (m=1:2:L-1) 
    signs(m)   = -1.0;
    signs(m+1) = 1.0;
end 
% Skip the s_00 component as it is zero 
for el = 1:L-1 
    % N.B. for curvelets, m = el; 
    % N.B. the condition : m < L is satisfied;
    m = el; 
    % for positive m
    ind_pm = ssht_elm2ind(el, m);
    s_lm(ind_pm)= sqrt(1./2.);
    % for negative m
    ind_nm = ssht_elm2ind(el, -m);
    s_lm(ind_nm)= signs(m)*conj(s_lm(ind_pm));
end


% ----------
% Curvelet angular components:
% ----------
[kappa kappa0] =  s2let_transform_axisym_tiling(B, L, J_min);
el_min = max(abs(Spin), abs(original_spin));
for j = J_min:J
 for el = el_min:L-1  
     m = el; 
     % positive m
     ind_pm = ssht_elm2ind(el, m); 
     cur_lm{j-J_min+1}(ind_pm) = s_lm(ind_pm) * sqrt((2*el+1)/(8.0*pi*pi))* kappa(j+1,el+1) ; 
     % negative m
     ind_nm = ssht_elm2ind(el, -m); 
     cur_lm{j-J_min+1}(ind_nm) =((-1)^m)* conj(cur_lm{j-J_min+1}(ind_pm)) ; 
     % if SpinLowered == true 
     if (args.SpinLowered ~= 0)  
      s2let_spin_lowered_norm_factor = s2let_spin_lowered_normalization(el, 'original_spin',original_spin);
      cur_lm{j-J_min+1}(ind_pm) = cur_lm{j-J_min+1}(ind_pm)*s2let_spin_lowered_norm_factor ;
      cur_lm{j-J_min+1}(ind_nm) = cur_lm{j-J_min+1}(ind_nm)*s2let_spin_lowered_norm_factor ;
     end 
 end
end

% ----------
% Scaling Function
% ----------
scal_l = zeros(L^2,1);
for el = el_min:L-1 
    scal_l(el^2+el+1,1) = sqrt((2*el+1)/(4.0*pi)) *kappa0(el+1); 
    % if SpinLowered == true    
    if (args.SpinLowered ~= 0)  
     s2let_spin_lowered_norm_factor = s2let_spin_lowered_normalization(el, 'original_spin',original_spin);
     scal_l(el^2+el+1,1) = scal_l(el^2+el+1,1) *s2let_spin_lowered_norm_factor ;
    end 
end


end
