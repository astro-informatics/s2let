function [f_lmn, f_scal_lm] = s2let_transform_analysis_lm2lmn(flm, varargin)

% s2let_transform_analysis_lm2cur
% Compute spin directional curvelet transform, input in harmonic space,
% output in pixel space.
%
% Default usage :
%
%   [f_cur, f_scal] = s2let_transform_analysis_lm2cur(flm, <options>)
%
% flm is the input field in harmonic space,
% f_cur contains the output curvelet contributions,
% f_scal contains the output scaling contributions,
%
% Option :
%  'B'               = { Dilation factor; B > 1 (default=2) }
%  'L'               = { Harmonic band-limit; L > 1 (default=guessed from input) }
%  'N'               = { Azimuthal/directional band-limit; N > 1 (default=L) }
%  'Spin'               = { Spin; (default=0) }
%  'J_min'           = { Minimum curvelet scale to consider;
%                        0 <= J_min < log_B(L) (default=0) }
%  'Upsample'      = { false        [multiresolution algorithm (default)],
%                      true       [full resolution curvelets] }
%  'Sampling'        = { 'MW'           [McEwen & Wiaux sampling (default)],
%                        'MWSS'         [McEwen & Wiaux symmetric sampling] }
%  'Reality'         = { false        [do not assume corresponding signal f real (default)],
%                        true         [assume f real (improves performance)] }
%  'SpinLowered'     = { true  [Apply normalisation factors for spin-lowered
%                               curvelets and scaling function.],
%                        false [Apply the usual normalisation factors such
%                               that the curvelets fulfil the admissibility
%                               condition (default)]}
%  'SpinLoweredFrom' = [integer; if the SpinLowered option is used, this
%                       option indicates which spin number the curvelets
%                       should be lowered from (default = 0)]
%
% S2LET package to perform curvelets transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

sz = length(flm(:));
Lguessed = sqrt(sz);

p = inputParser;
p.addRequired('flm', @isnumeric);
p.addParamValue('B', 2, @isnumeric);
p.addParamValue('L', Lguessed, @isnumeric);
p.addParamValue('J_min', 0, @isnumeric);
p.addParamValue('N', Lguessed, @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('Upsample', false, @islogical);
p.addParamValue('Sampling', 'MW', @ischar);
p.addParamValue('Reality', false, @islogical);
p.addParamValue('SpinLowered', false, @islogical);
p.addParamValue('SpinLoweredFrom', 0, @isnumeric);
p.parse(flm, varargin{:});
args = p.Results;

flm_vec = flm(:);

if(all(isreal(flm_vec)))
  flm_vec = complex(flm_vec,0);
end

[f_cur_vec, f_scal_vec] = s2let_transform_analysis_lm2cur_mex(flm_vec, args.B, args.L, args.J_min, ...
                                                              args.N, args.Spin, ...
                                                              args.Reality, args.Upsample, ...
                                                              args.SpinLowered, args.SpinLoweredFrom, ...
                                                              args.Sampling);


% Generate flmn of complex signal (c.f. s2let_analysis_lm2lmn(f_cur_lmn, f_scal_lm, flm, cur_lm, scal_l, parameters))

ind_ln=0;
ind = 0;
ind_lmn = 0;
for j = J_min:J,
flmn{j-J_min+1} = zeros((2*N-1)*L*L,1);
for n = -N+1:N-1,
for el = abs(n):L-1,
ind_ln = ssht_elm2ind(el, n);   %for l and n
psi = 8.*pi*pi/(2.*el+1) *conj(cur_lm{j-J_min+1}(ind_ln));
for m = -el:el,
ind = ssht_elm2ind(el, m);
ind_lmn = so3_elmn2ind(el,m,n,L,N);
flmn{j-J_min+1}(ind_lmn) =  flm_gen(ind) * psi;
end
end
end
end
%% Scaling function
disp('Compute scaling function f_scal_lm=flm_gen(lm_ind) * phi ')
%f_scal_lm = zeros(L^2,1)
lm_ind=0;
for j = J_min:J,
for el = 0:L-1,
phi = sqrt(4.0*pi/(2.*el+1)) * phi_l(el+1);  %
for m = -el:el,
lm_ind=ssht_elm2ind(el, m);
f_scal_lm(lm_ind) = flm_gen(lm_ind) * phi;
end
end
end

% Compute inverse then forward transform.
disp('so3_forward (i.e. f to flmn) and so3_inverse (i.e. flmn to f): ')
for j = J_min:J,
f_cur = so3_inverse(flmn{j-J_min+1}, L, N);       %i.e. so3_core_inverse_via_ssht -> f_cur[nalpha*nbeta*ngamma*sizeof(f)]
flmn_syn{j-J_min+1} = so3_forward(f_cur, L, N);   %i.e. so3_core_forward_via_ssht -> input funcion on the sphere and output flmn harmonic coefficient
% Compute maximum error in harmonic space.
disp('Check the error of so3_forward (i.e. f to flmn) and so3_inverse (i.e. flmn to f): ')
maxerr = max(abs(flmn_syn{j-J_min+1} - flmn{j-J_min+1}))
end