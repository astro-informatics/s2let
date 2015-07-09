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
J = s2let_jmax(args.L, args.B);

% ------------
% Signal Analysis: 
% ------------
% Generate flmn from flm of the complex signal
% (c.f. C function: s2let_analysis_lm2lmn(f_cur_lmn, f_scal_lm, flm, cur_lm, scal_l, parameters);
% -----------------
% Curvelet contribution:
% -----------------
% disp('ana_lm2lmn: Curvelet analysis of complex signals from harmonic to Wigner space (i.e. flm to flmn)')
 
for j = args.J_min:J,   
  band_limit = min([ s2let_bandlimit(j,args.J_min,args.B,args.L) args.L ]);
  Nj = band_limit;
  if (args.Upsample ~= 0)  %true => full-resolution, set band_limit back to args.L
      band_limit = args.L;
  end
% for the case SO3_STORAGE_PADDED:
  if (args.Reality == 0) %i.e. false (default) => complex signals
      f_cur_lmn{j-args.J_min+1} = zeros((2*Nj-1)*band_limit^2,1);
  else %i.e. real signals
      f_cur_lmn{j-args.J_min+1} = zeros(Nj*band_limit^2,1);  
  end 
%{
% TODO: in the future add option for SO3_STORAGE: 
% for the case SO3_STORAGE_COMPACT:
  if ((args.Reality == 0) && (args.Storage =='Compact') )
  f_cur_lmn{j-args.J_min+1} = zeros((2*Nj-1)*(3*band_limit^2-Nj*(Nj-1)/3,1);
  else if ((args.Reality == 1) && (args.Storage =='Compact') )
  f_cur_lmn{j-args.J_min+1} = zeros(Nj*(6*band_limit^2-(Nj-1)*(2*Nj-1))/6,1);  
  end 
%} 
  ind_ln=0;
  ind_lm = 0;
  ind_lmn = 0;
  %
  if (args.Reality == 0) %i.e. false (default) => complex
   for en = -Nj+1:Nj-1,
       for el = max(abs(args.Spin),abs(en)):band_limit-1,
           ind_ln = ssht_elm2ind(el, en);
           psi = 8.*pi*pi/(2.*el+1) *conj(cur_lm{j-args.J_min+1}(ind_ln));
           for m = -el:el,
               ind_lm = ssht_elm2ind(el, m);
               if (args.Upsample == 0)  %false => multi-resolution
                        ind_lmn = so3_elmn2ind(el,m,en,band_limit,Nj);
               else
                        ind_lmn = so3_elmn2ind(el,m,en,args.L,Nj);
               end % end the if-loop for upsample
               f_cur_lmn{j-args.J_min+1}(ind_lmn) =  flm_init(ind_lm) * psi;
           end
       end
   end
  else % i.e.(args.Reality == 1) %i.e. true => real
   for en = 1-mod(Nj,2):Nj-1,
       for el = en:band_limit-1,
           ind_ln = ssht_elm2ind(el, en);
           psi = 8.*pi*pi/(2.*el+1) *conj(cur_lm{j-args.J_min+1}(ind_ln));
           for m = -el:el,     
               ind_lm = ssht_elm2ind(el, m);
               ind_lmn = so3_elmn2ind(el,m,en,band_limit,Nj,'Reality', args.Reality);
               f_cur_lmn{j-args.J_min+1}(ind_lmn) =  flm_init(ind_lm) * psi;
           end
       end
   end
  end % end if loop for Reality Option
end % end j-loop 
% For debugging: 
% disp('--check the size of f_cur_lmn--')
 % (2N-1)*L^2 for complex signals -> so3_storage_padded
 %  N*L^2 for real signals  -> so3_storage_padded
%whos f_cur_lmn
%len= length(f_cur_lmn)
%temp = f_cur_lmn{len};
%sz = size(temp)
%disp('--------')
% For L=N=16, real data: sz= 4096 1: 
   
% -----------------
% Scaling function contribution: 
% -----------------
% disp('ana_lm2lmn: : Compute scaling function f_scal_lm=flm_init(lm_ind) * phi '); 
if (args.Upsample == 0)   %i.e. false (default) => multi-resolution 
    band_limit = min([ s2let_bandlimit(args.J_min-1, args.J_min, args.B,args.L) args.L ]);
else 
    band_limit = args.L ; 
end
f_scal_lm = zeros(band_limit^2,1);
% ToDO Check: band_limit*(2*band_limit-1) for MW; (band_limit+1)*2*band_limit for MWSS
lm_ind=0;
if (args.Reality == 0) %i.e. false (default) => complex
   for el = abs(args.Spin):band_limit-1,
       phi = sqrt(4.0*pi/(2.*el+1))*scal_l(el^2+el+1,1);
       for m = -el:el,
           lm_ind=ssht_elm2ind(el, m);
           f_scal_lm(lm_ind) = flm_init(lm_ind) * phi;
       end
   end
else   %i.e. (args.Reality == 1); true => real
   for el = 0 :band_limit-1,
       phi = sqrt(4.0*pi/(2.*el+1))*scal_l(el^2+el+1,1);
       for m = -el:el,
           lm_ind=ssht_elm2ind(el, m);
           f_scal_lm(lm_ind) = flm_init(lm_ind) * phi;
       end
   end
end  % end if-loop for Reality Option
   
   
   
end

