function [f_cur_lmn, f_scal_lm] = s2let_transform_analysis_lm2lmn(flm_gen, cur_lm, scal_l, varargin)

% s2let_transform_analysis_lm2lmn
% Compute (spin) curvelet transform, input in harmonic space,
% output in Wigner space.
%
% Default usage :
%
%   [f_cur_lmn, f_scal_lm] = s2let_transform_analysis_lm2lmn(flm, cur_lm, scal_l, <options>)
%
% flm is the input field in harmonic space,
% cur_lm contains the curvelet kernels in harmonic space,
% scal_l contains the scaling contributions in harmonic space,
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

sz = length(flm_gen(:));
Lguessed = sqrt(sz);

p = inputParser;
p.addRequired('flm_gen', @isnumeric);
p.addRequired('cur_lm', @iscell);
p.addRequired('scal_l', @isnumeric);
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
p.parse(flm_gen,  cur_lm, scal_l, varargin{:});
args = p.Results;


%
% 'Signal Analysis: '
% Generate flmn from flm of the complex signal
% (c.f. C function: s2let_analysis_lm2lmn(f_cur_lmn, f_scal_lm, flm, cur_lm, scal_l, parameters);
%
disp('Curvelet analysis of complex signals in Wigner space (i.e. flm to flmn)')
ind_ln=0;
ind = 0;
ind_lmn = 0;
J = s2let_jmax(args.L, args.B);
for j = args.J_min:J,
  f_cur_lmn{j-args.J_min+1} = zeros((2*args.N-1)*args.L*args.L,1);
  for n = -args.N+1:args.N-1,
   for el = abs(n):args.L-1,
     ind_ln = ssht_elm2ind(el, n);
     psi = 8.*pi*pi/(2.*el+1) *conj(cur_lm{j-args.J_min+1}(ind_ln));
     for m = -el:el,
      ind = ssht_elm2ind(el, m);
      ind_lmn = so3_elmn2ind(el,m,n,args.L,args.N);
      f_cur_lmn{j-args.J_min+1}(ind_lmn) =  flm_gen(ind) * psi;
     end
    end
  end
end


%
% Scaling function
%
disp('Compute scaling function f_scal_lm=flm_gen(lm_ind) * phi '); 
f_scal_lm = zeros(args.L^2,1);
lm_ind=0;
for j = args.J_min:J,
 for el = 0:args.L-1,
   phi = sqrt(4.0*pi/(2.*el+1))*scal_l(el^2+el+1,1);     % phi_l(l+1);
   for m = -el:el,
    lm_ind=ssht_elm2ind(el, m);
    f_scal_lm(lm_ind) = flm_gen(lm_ind) * phi;
   end
 end
end

