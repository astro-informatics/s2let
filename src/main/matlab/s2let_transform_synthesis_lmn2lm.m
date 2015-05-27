function [flm_cur_syn, flm_scal_syn] = s2let_transform_synthesis_lmn2lm(f_cur_lmn, f_scal_lm, cur_lm, scal_l, flm_gen, varargin)

% s2let_transform_synthesis_lmn2lm
% Compute (spin) curvelet transform, input in Wigner space,
% output in harmonic space.
%
% Default usage :
%
%   [flm_cur_syn, f_scal_syn] = s2let_transform_synthesis_lmn2lm(f_cur_lmn, f_scal_lm, cur_lm, scal_l, <options>)
%
% f_cur_lmn is the input curvelet contribution,
% f_cur_lmn is the input scaling function contribution ,
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
p.addRequired('f_cur_lmn', @iscell);
p.addRequired('f_scal_lm', @isnumeric);
p.addRequired('cur_lm', @iscell);
p.addRequired('scal_l', @isnumeric);
p.addRequired('flm_gen', @isnumeric);
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
p.parse(f_cur_lmn, f_scal_lm, cur_lm, scal_l, flm_gen, varargin{:});
args = p.Results;


%
% 'Signal synthesis: '
% Generate flmn from flm of the complex signal
% (c.f. C function: s2let_synthesis_cur_lmn2lm(flm, f_cur_lmn, f_scal_lm, cur_lm, scal_l, parameters);)
%
disp('Compute flm_cur_syn from flmn_syn');
ind_ln =0;
ind=0;
ind_lmn=0;
flm_cur_syn=zeros(args.L^2,1);
J = s2let_jmax(args.L, args.B);
for j = args.J_min:J,
 for n = -args.N+1:args.N-1,
  for el = abs(n):args.L-1,
   ind_ln = ssht_elm2ind(el, n);
   psi = (cur_lm{j-args.J_min+1}(ind_ln));
   for m = -el:el,
    ind = ssht_elm2ind(el, m);
    ind_lmn = so3_elmn2ind(el,m,n,args.L,args.N);
    flm_cur_syn(ind) =flm_cur_syn(ind)+ f_cur_lmn{j-args.J_min+1}(ind_lmn)* psi;
   end
  end
 end
end

%
% Scaling function
%
disp('Compute flm_scal_syn ');
flm_scal_syn=zeros(args.L^2,1);
lm_ind=0;
for el = 0:args.L-1,
  phi = sqrt(4.*pi/(2.*el+1))*scal_l(el^2+el+1,1);         % phi_l(l+1);
  for m = -el:el,
   lm_ind=ssht_elm2ind(el, m);
   flm_scal_syn(lm_ind) =  flm_scal_syn(lm_ind)+ flm_gen(lm_ind)* phi;
  end
end







