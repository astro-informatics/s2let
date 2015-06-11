function flm_rec= s2let_transform_curvelet_synthesis_lmn2lm(f_cur_lmn, f_scal_lm, cur_lm,  scal_l, varargin)

% s2let_transform_synthesis_lmn2lm
% Compute (spin) curvelet transform, 
% input in Wigner space,
% output in harmonic space.
%
% Default usage :
%
%   flm_rec = s2let_transform_curvelet_synthesis_lmn2lm(f_cur_lmn, f_scal_lm,  <options>)
%
% f_cur_lmn is the input curvelet contribution,
% f_cur_lmn is the input scaling function contribution,
% cur_lm is the curvelet kernels, 
% scal_l is the scaling function kernel. 
%
% Option :
%  'B'               = { Dilation factor; B > 1 (default=2) }
%  'L'               = { Harmonic band-limit; L > 1 (default=guessed from input) }
%  'Spin'               = { Spin; (default=0) }
%  'J_min'           = { Minimum curvelet scale to consider;
%                        0 <= J_min < log_B(L) (default=0) }
%  'Upsample'      = { false        [multiresolution algorithm (default)],
%                      true       [full resolution curvelets] }
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
%  'Sampling'        = { 'MW'           [McEwen & Wiaux sampling (default)],
%                        'MWSS'         [McEwen & Wiaux symmetric sampling] }
%
% -----------------------------------------------------------
% Log: 
% -  constructed by Jennifer Y H Chan on 5th June 2015  
% -----------------------------------------------------------
% S2LET package to perform wavelet transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details
% -----------------------------------------------------------
p = inputParser;
p.addRequired('f_cur_lmn', @iscell);
p.addRequired('f_scal_lm', @isnumeric);
p.addRequired('cur_lm', @iscell);
p.addRequired('scal_l', @isnumeric);
p.addParamValue('B', 2, @isnumeric);
p.addParamValue('L',  @isnumeric);   
p.addParamValue('J_min', 0, @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('Upsample', false, @islogical);
p.addParamValue('Sampling', 'MW', @ischar);
p.addParamValue('Reality', false, @islogical);
p.addParamValue('SpinLowered', false, @islogical);
p.addParamValue('SpinLoweredFrom', 0, @isnumeric);
p.parse(f_cur_lmn, f_scal_lm, cur_lm, scal_l, varargin{:});
args = p.Results;

N = args.L ;
Nj = N; 
band_limit = args.L; 

% Debug:
% f_scal_plot = ssht_inverse(f_scal_lm, args.L, 'Reality', args.Reality);
% ssht_plot_mollweide(f_scal_plot, args.L, 'Mode', 1);

% ---------------
% Signal synthesis: 
% ---------------
% Generate flmn from flm of the complex signal
% (c.f. C function: s2let_synthesis_cur_lmn2lm(flm, f_cur_lmn, f_scal_lm, cur_lm, scal_l, parameters);)
% -----------------
% Curvelet contribution:
% -----------------
% disp('syn_lmn2lm: Curvelet synthesis of complex signals from Wigner to harmonic space (flm_cur_syn from flmn_syn)');

flm_cur_syn=zeros(args.L^2,1);
J = s2let_jmax(args.L, args.B);
for j = args.J_min:J,
 ind_ln =0;
 ind=0;
 ind_lmn=0;
 if (args.Upsample ~= 1)  %false => multi-resolution 
     band_limit = min([ s2let_bandlimit(j,args.J_min,args.B,args.L) args.L ]);
     Nj =  band_limit;
     %Nj = min(N, band_limit);
     %Nj = Nj+ mod((Nj+N),2) ;  %ensure N and Nj are both even and both odd
 end
 for en = -Nj+1:Nj-1,
  for el = max(abs(args.Spin),abs(en)):band_limit-1,
   ind_ln = ssht_elm2ind(el, en);
   psi = (cur_lm{j-args.J_min+1}(ind_ln));
   for m = -el:el,
    ind = ssht_elm2ind(el, m);
    ind_lmn = so3_elmn2ind(el,m,en,band_limit,Nj);
    flm_cur_syn(ind) =flm_cur_syn(ind)+ f_cur_lmn{j-args.J_min+1}(ind_lmn)* psi;
   end
  end
 end
end

% -----------------
% Scaling function contribution: 
% -----------------
% disp('syn_lmn2lm:  Compute flm_scal_syn ');
if (args.Upsample ~= 1)  %false => multi-resolution 
   band_limit = min([ s2let_bandlimit(args.J_min-1, args.J_min, args.B,args.L) args.L ]);
end
flm_scal_syn=zeros(args.L^2,1);
lm_ind=0;
for el = abs(args.Spin): band_limit-1,
  phi = sqrt(4.*pi/(2.*el+1))*scal_l(el^2+el+1,1);      
  for m = -el:el,
   lm_ind=ssht_elm2ind(el, m);
   flm_scal_syn(lm_ind) =  flm_scal_syn(lm_ind)+  f_scal_lm(lm_ind)* phi;
  end
end

% -----------------
% Summing over their contributions: 
% -----------------
% disp('Sum: flm_cur_syn+flm_scal_syn ');
flm_rec = flm_scal_syn + flm_cur_syn;

% Clear array memory: 
flm_init =0;

end






