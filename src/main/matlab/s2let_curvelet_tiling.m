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
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

p = inputParser;
p.addRequired('B', @isnumeric);
p.addRequired('J_min', @isnumeric);
p.addRequired('L',  @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('SpinLowered', false, @islogical);
p.addParamValue('SpinLoweredFrom', 0, @isnumeric);
p.parse(B, L, J_min, varargin{:});
args = p.Results;

J = s2let_jmax(L, B);

[kappa kappa0] =  s2let_transform_axisym_tiling(B, L, J_min);

% Curvelet Kernels:
ind_pm = 0;
ind_nm = 0;
for j = J_min:J
 for el = 0:L-1
 m = el;
 % for positive m
 ind_pm = ssht_elm2ind(el, m);
 % Curvelet coefficients:
 cur_lm{j-J_min+1}(ind_pm) = kappa(j+1,el+1);
 % for negative m
 ind_nm = ssht_elm2ind(el, -m);
 cur_lm{j-J_min+1}(ind_nm) = (-1)^m * conj(cur_lm{j-J_min+1}(ind_pm));
  end
end 

% Scaling Function
scal_l = zeros(L^2,1);
for l = 0:L-1
 scal_l(l^2+l+1,1) = kappa0(l+1);
end



end
