function [f_cur_lmn, f_scal_lm] = s2let_transform_curvelet_analysis_lm2lmn(flm_init, cur_lm, scal_l, varargin)

% s2let_transform_curvelet_analysis_lm2lmn
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
%  'Spin'            = { Spin; (default=0) }
%  'Reality'         = { false        [do not assume corresponding signal f real (default)],
%                        true         [assume f real (improves performance)] }
%  'Upsample'        = { false        [multiresolution algorithm (default)],
%                        true         [full resolution curvelets] }
%  'SpinLowered'     = { true  [Apply normalisation factors for spin-lowered
%                               curvelets and scaling function.],
%                        false [Apply the usual normalisation factors such
%                               that the curvelets fulfil the admissibility
%                               condition (default)]}
%  'SpinLoweredFrom' = [integer; if the SpinLowered option is used, this
%                       option indicates which spin number the curvelets
%                       should be lowered from (default = 0)]
% -----------------------------------------------------------
% S2LET package to perform Wavelet Transform on the Sphere.
% Copyright (C) 2012-2016  Boris Leistedt, Jennifer Chan & Jason McEwen
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

J = s2let_jmax(args.L, args.B);

% -----------------
% Signal Analysis: 
% -----------------
% Curvelet contribution:  
for j = args.J_min:J,   
  band_limit = min([ s2let_bandlimit(j,args.J_min,args.B,args.L) args.L ]);
  % Nj = orientational band-limit at j-th scale 
  Nj = band_limit;
  if (args.Reality == 0) 
      % for the case SO3_STORAGE_PADDED:
      if (args.Upsample ~= 0) 
           f_cur_lmn{j-args.J_min+1} = zeros((2*Nj-1)*args.L^2,1); 
      else 
           f_cur_lmn{j-args.J_min+1} = zeros((2*Nj-1)*band_limit^2,1);
      end 
      for en = -Nj+1:Nj-1,
          for el = max(abs(args.Spin),abs(en)):band_limit-1,
              ind_ln = ssht_elm2ind(el, en);
              psi = 8.*pi*pi/(2.*el+1) *conj(cur_lm{j-args.J_min+1}(ind_ln));
              for m = -el:el,
                  ind_lm = ssht_elm2ind(el, m);
                  if (args.Upsample ~= 0)  
                      ind_lmn = so3_elmn2ind(el,m,en,args.L,Nj);
                  else
                      ind_lmn = so3_elmn2ind(el,m,en,band_limit,Nj);
                  end 
                  f_cur_lmn{j-args.J_min+1}(ind_lmn) =  flm_init(ind_lm) * psi;
               end
          end
      end
  else % real setting: 
      if (args.Upsample ~= 0) 
          f_cur_lmn{j-args.J_min+1} = zeros(Nj*args.L^2,1);  
      else 
          f_cur_lmn{j-args.J_min+1} = zeros(Nj*band_limit^2,1);  
      end 
      for el = 1:band_limit-1,
          for m = -el:el, 
              ind_lm = ssht_elm2ind(el, m);
              % (n =0) terms:  
              ind_lnzero = ssht_elm2ind(el, 0);
              psizero = 8.*pi*pi/(2.*el+1) *conj(cur_lm{j-args.J_min+1}(ind_lnzero));
              ind_lmnzero = so3_elmn2ind(el,m,0,band_limit,Nj,'Reality', args.Reality);
              f_cur_lmn{j-args.J_min+1}(ind_lmnzero) =  flm_init(ind_lm) *psizero; 
              % (n ~=0) terms:  
              for en = 1: Nj-1,    
                  if (el >= en)
                      ind_ln = ssht_elm2ind(el, en);
                      psi = 8.*pi*pi/(2.*el+1) *conj(cur_lm{j-args.J_min+1}(ind_ln));
                      if (args.Upsample == 0) 
                          ind_lmn = so3_elmn2ind(el,m,en,band_limit,Nj,'Reality', args.Reality);
                      else
                          ind_lmn = so3_elmn2ind(el,m,en,args.L,Nj,'Reality', args.Reality);
                      end 
                      f_cur_lmn{j-args.J_min+1}(ind_lmn) =  flm_init(ind_lm) * psi;
                  end
              end 
          end
      end
  end % end if loop for Reality Option
end % end j-loop 


% Scaling function contribution: 
if (args.Upsample ~= 0)  
    band_limit = args.L ; 
else 
    band_limit = min([ s2let_bandlimit(args.J_min-1, args.J_min, args.B,args.L) args.L ]);
end
f_scal_lm = zeros(band_limit^2,1);
if (args.Reality == 0)
   for el = abs(args.Spin):band_limit-1,
       phi = sqrt(4.0*pi/(2.*el+1))*scal_l(el^2+el+1,1);
       for m = -el:el,
           lm_ind=ssht_elm2ind(el, m);
           f_scal_lm(lm_ind) = flm_init(lm_ind) * phi;
       end
   end
else  % real setting: 
   for el = 0 :band_limit-1,
       phi = sqrt(4.0*pi/(2.*el+1))*scal_l(el^2+el+1,1);
       for m = -el:el,
           lm_ind=ssht_elm2ind(el, m);
           f_scal_lm(lm_ind) = flm_init(lm_ind) * phi;
       end
   end
end  % end if-loop for Reality Option
   
   
end

