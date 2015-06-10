function [f_cur_lmn, f_scal_lm] = s2let_transform_curvelet_analysis_lm2lmn(flm_init, cur_lm, scal_l, varargin)

% s2let_transform_analysis_lm2lmn
% Compute (spin) curvelet transform, input in harmonic space,
% output in Wigner space.
%
% Default usage :
%
%   [f_cur_lmn, f_scal_lm] = s2let_transform_curvelet_analysis_lm2lmn(flm_init, cur_lm, scal_l, <options>)
%
% flm_init is the input field in harmonic space,
% cur_lm is the curvelet kernels, 
% scal_l is the scaling function kernel. 
%
% Option :
%  'B'               = { Dilation factor; B > 1 (default=2) }
%  'L'               = { Harmonic band-limit; L > 1 (default=guessed from input) }
%  'J_min'           = { Minimum curvelet scale to consider;
%                        0 <= J_min < log_B(L) (default=0) }
%  'Spin'               = { Spin; (default=0) }
%  'Reality'         = { false        [do not assume corresponding signal f real (default)],
%                        true         [assume f real (improves performance)] }
%  'Upsample'      = { false        [multiresolution algorithm (default)],
%                      true         [full resolution curvelets] }
%  'SpinLowered'     = { true  [Apply normalisation factors for spin-lowered
%                               curvelets and scaling function.],
%                        false [Apply the usual normalisation factors such
%                               that the curvelets fulfil the admissibility
%                               condition (default)]}
%  'SpinLoweredFrom' = [integer; if the SpinLowered option is used, this
%                       option indicates which spin number the curvelets
%                       should be lowered from (default = 0)]
%  'Sampling'        = { 'MW'           [McEwen & Wiaux sampling (default)],
%                        'MWSS'         [McEwen & Wiaux symmetric sampling] }
% -----------------------------------------------------------
% Log: 
% -  constructed by Jennifer Y H Chan on 5th June 2015  
% -----------------------------------------------------------
% S2LET package to perform wavelet transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details
% -----------------------------------------------------------

sz = length(flm_init(:));
Lguessed = sqrt(sz);

p = inputParser;
p.addRequired('flm_init', @isnumeric);
p.addRequired('cur_lm', @iscell);
p.addRequired('scal_l', @isnumeric);
p.addParamValue('B', 2, @isnumeric);
p.addParamValue('L', Lguessed, @isnumeric);
p.addParamValue('J_min', 0, @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('Reality', false, @islogical);
p.addParamValue('Upsample', false, @islogical);
p.addParamValue('SpinLowered', false, @islogical);
p.addParamValue('SpinLoweredFrom', 0, @isnumeric);
p.addParamValue('Sampling', 'MW', @ischar);
p.parse(flm_init, cur_lm, scal_l, varargin{:});
args = p.Results;

% Azimuthal/directional band-limit N =L always holds for curvelets:           
N = args.L ;
% Nj = N; 
% band_limit = args.L; 
J = s2let_jmax(args.L, args.B);

% ---------------
% Tile curvelets:
% ---------------
%[cur_lm scal_l] = s2let_curvelet_tiling(args.B, args.L, args.J_min, ...
%                                        'Spin', args.Spin, 'SpinLowered', args.SpinLowered,...
%                                        'SpinLoweredFrom',args.SpinLoweredFrom);
% where 
% cur_lm, @iscell, contains the curvelet kernels in harmonic space,
% scal_l, @isnumeric, contains the scaling contributions in harmonic space


% ------------
% Signal Analysis: 
% ------------
% Generate flmn from flm of the complex signal
% (c.f. C function: s2let_analysis_lm2lmn(f_cur_lmn, f_scal_lm, flm, cur_lm, scal_l, parameters);
% -----------------
% Curvelet contribution:
% -----------------
% disp('ana_lm2lmn: Curvelet analysis of complex signals from harmonic to Wigner space (i.e. flm to flmn)')

%{
  if (args.Upsample == 0)  %i.e. false => multi-resolution
    band_limit = min([ s2let_bandlimit(j,args.J_min,args.B,args.L) args.L ]);
    %Nj = min(N, band_limit);
    %Nj = Nj+ mod((Nj+N),2) ;  % ensure Nj and J are both even or both odd
    Nj = band_limit;
  end
%}  
ind_ln=0;
ind = 0;
ind_lmn = 0;
for j = args.J_min:J,
  f_cur_lmn{j-args.J_min+1} = zeros((2*N-1)*args.L*args.L,1);
  for en =  -N+1: N+1,   
       %-Nj+1:Nj-1,
   for el = max(abs(args.Spin),abs(en)):args.L-1,   
       %band_limit-1,
     ind_ln = ssht_elm2ind(el, en);
     psi = 8.*pi*pi/(2.*el+1) *conj(cur_lm{j-args.J_min+1}(ind_ln));
     for m = -el:el,
      ind = ssht_elm2ind(el, m);
      ind_lmn = so3_elmn2ind(el,m,en,args.L,N);   %band_limit,Nj);
      f_cur_lmn{j-args.J_min+1}(ind_lmn) =  flm_init(ind) * psi;
     end
    end
  end
end


% -----------------
% Scaling function contribution: 
% -----------------
% disp('ana_lm2lmn: : Compute scaling function f_scal_lm=flm_init(lm_ind) * phi '); 
%{
if (args.Upsample ~= 1)  %false => multi-resolution 
   band_limit = min([ s2let_bandlimit(args.J_min-1, args.J_min, args.B,args.L) args.L ]);
end
%}
f_scal_lm = zeros(args.L^2,1);
lm_ind=0;
for el =  0:args.L-1, %abs(args.Spin):band_limit-1, 
   phi = sqrt(4.0*pi/(2.*el+1))*scal_l(el^2+el+1,1);   
   for m = -el:el,
    lm_ind=ssht_elm2ind(el, m);
    f_scal_lm(lm_ind) = flm_init(lm_ind) * phi;
   end
end


% Clear array memory: 
flm_init =0;

end

