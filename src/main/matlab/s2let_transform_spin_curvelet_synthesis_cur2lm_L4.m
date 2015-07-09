function flm_rec = s2let_transform_spin_curvelet_synthesis_cur2lm(f_cur, f_scal,  varargin)

% s2let_transform_spin_curvelet_synthesis_cur2lm
% Compute curvelet transform:  
% input in curvelet space (i.e. hamonic to Wigner space via SO3_forward) 
% output in harmonic space (i.e. Wigner to harmonic space via synthesis_lmn2lm) .
%
% Default usage :
%
% flm_rec = s2let_transform_spin_curvelet_synthesis_lm2cur(f_cur, f_scal,  <options>)
%
% f_cur contains the input curvelet contributions -- MW sampling,
% f_scal contains the input scaling contributions -- MW sampling,
% flm_rec is the output field = flm_cur_syn+ flm_scal_syn
%
% Option :
%  'B'               = { Dilation factor; B > 1 (default=2) }
%  'L'               = { Harmonic band-limit; L > 1 (default=guessed from input) }
%  'J_min'           = { Minimum curvelet scale to consider;
%  'Spin'               = { Spin; (default=0) }
%                        0 <= J_min < log_B(L) (default=0) }
%  'Reality'         = { false        [do not assume corresponding signal f real (default)],
%                        true         [assume f real (improves performance)] }
%  'Upsample'        = { false        [multiresolution algorithm (default)],
%                        true       [full resolution curvelets] }
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

len = size(f_cur);
temp = f_cur{len};
sz = size(temp);
Lguessed = sz(2);
%if sz(1) == 2*sz(2)-1 || sz(2) == 2*sz(1)-1
%    Lguessed = min([sz(1) sz(2)]) 
%else
%    Lguessed = min([sz(1) sz(2)])-1 ;
%end

p = inputParser;
p.addRequired('f_cur');
p.addRequired('f_scal', @isnumeric);
p.addParamValue('B', 2, @isnumeric);
p.addParamValue('L', Lguessed, @isnumeric);                           
p.addParamValue('J_min', 0, @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('Upsample', false, @islogical);
p.addParamValue('Reality', false, @islogical);
p.addParamValue('SpinLowered', false, @islogical);
p.addParamValue('SpinLoweredFrom', 0, @isnumeric);
p.addParamValue('Sampling', 'MW', @ischar);
p.parse(f_cur, f_scal, varargin{:});
args = p.Results;


N= args.L;
J = s2let_jmax(args.L, args.B);

% ---------------
% Tile curvelets:
% ---------------
% Call curvelet- and scaling-function- generating functions
[cur_lm scal_l] = s2let_spin_curvelet_tiling(args.B, args.L, args.J_min,...
                                            'Spin', args.Spin,...
                                            'SpinLowered', args.SpinLowered, ...
                                            'SpinLoweredFrom', args.SpinLoweredFrom);

% -----------------
% Signal synthesis: (Transform to lmn space, then reconstruct the signal in harmonic space)
% -----------------
% Scaling functions:
if (args.Upsample == 0)  %false => multi-resolution 
     band_limit = min([ s2let_bandlimit(args.J_min-1,args.J_min,args.B,args.L) args.L ]);
else
     band_limit = args.L ;
end
f_scal_lm_syn = ssht_forward(f_scal, band_limit, ...
                            'Method', args.Sampling,...
                            'Spin', 0, ...
                            'Reality', args.Reality);
  
% Curvelet contribution:
for j = args.J_min:J,
    band_limit = min([ s2let_bandlimit(j,args.J_min,args.B,args.L) args.L ]);
    Nj = band_limit;
    if (args.Upsample ==0)  %false => multi-resolution
        f_cur_lmn_syn{j-args.J_min+1} = so3_forward(f_cur{j-args.J_min+1} , band_limit, Nj, ...
                                                    'Reality', args.Reality, 'Sampling', args.Sampling);
    else  % Upsample true => full-resolution
        f_cur_lmn_syn{j-args.J_min+1} = so3_forward(f_cur{j-args.J_min+1} , args.L, Nj, ...
                                                    'Reality', args.Reality, 'Sampling', args.Sampling);
    end
end

% Rotate the Wigner coefficients f_cur_lmn (such the curvelets centred at the North pole)
% Exploit the property of curvelets that cur_ln = (cur_ll)*(delta_ln)
% such that cur_lmk_rotated = cur_lml*conj(Dlkl(evaluated at the desired rotation angle))
% ---------------
% Define Euler angles
% (for rotating the curvelets to the north pole):
% ---------------
alpha = pi ;
beta = pi/2 ;
gamma = 0 ;
% ---------------
% Precompute Wigner small-d functions
% denoted here as d (in the paper: d_lmn for all el, m, n evaluated at beta).
% They are indexed d(el,m,n).
% Alpha and gamma are the other two rotation angles.
% ---------------
d = zeros(args.L, 2*args.L-1, 2*args.L-1);
d(1,:,:) = ssht_dl(squeeze(d(1,:,:)), args.L, 0, beta);
for el = 1:args.L-1
d(el+1,:,:) = ssht_dl(squeeze(d(el,:,:)), args.L, el, beta);
end

for j = args.J_min:J,
  band_limit = min([ s2let_bandlimit(j,args.J_min,args.B,args.L) args.L ]);
  Nj = band_limit; 
  % for the case SO3_STORAGE_PADDED:
  if (args.Reality == 0) %i.e. false (default) => complex signals
      if (args.Upsample ~= 0) 
          f_cur_lmn_syn_rotated{j-args.J_min+1} = zeros((2*Nj-1)*args.L^2,1); 
      else 
          f_cur_lmn_syn_rotated{j-args.J_min+1} = zeros((2*Nj-1)*band_limit^2,1);
      end 
  else %i.e. real signals
      if (args.Upsample ~= 0) 
          f_cur_lmn_syn_rotated{j-args.J_min+1} = zeros(Nj*args.L^2,1);  
      else
          f_cur_lmn_syn_rotated{j-args.J_min+1} = zeros(Nj*band_limit^2,1);  
      end
  end 
  ind_lmk = 0;
  ind_lmn = 0;
  for el = abs(args.Spin):(band_limit-1) %0: args.L-1  % max(abs(args.Spin),abs(en))
      for m = -el:el
          en_max = min(el, Nj-1); 
          if (args.Reality == 0) %i.e. false (default) => complex
             for k = -en_max:en_max
                 if (args.Upsample == 0)  %false => multi-resolution
                         ind_lmk = so3_elmn2ind(el,m,k,band_limit,Nj);
                 else
                         ind_lmk = so3_elmn2ind(el,m,k,args.L,Nj);
                 end % end the if-loop for upsample
                 for en = -en_max:en_max
                     %  Dlmn = exp(-1i*m*alpha) * d(el+1,m+L,n+L) * exp(-1i*n*gamma);
                     Dlkn = exp(-1i*k*alpha) * d(el+1,k+args.L,en+args.L) * exp(-1i*en*gamma);
                     if (args.Upsample == 0)  %false => multi-resolution
                         ind_lmn = so3_elmn2ind(el,m,en,band_limit,Nj);
                     else
                         ind_lmn = so3_elmn2ind(el,m,en,args.L,Nj);
                     end % end the if-loop for upsample
                     f_cur_lmn_syn_rotated{j-args.J_min+1}(ind_lmk)=  f_cur_lmn_syn_rotated{j-args.J_min+1}(ind_lmk)+ ...
                                                                      conj(Dlkn) * f_cur_lmn_syn{j-args.J_min+1}(ind_lmn);
                 end % end n-loop 
             end % end k-loop
          end % end if-loop for reality option 
      end % end m-loop 
  end % end el-loop                    
end

% Reconstruct the signals in harmonic space
% Perform Wigner transform (lmn2lm -  Call matlab function synthesis_lmn2lm)
flm_rec = s2let_transform_spin_curvelet_synthesis_lmn2lm(f_cur_lmn_syn_rotated, f_scal_lm_syn, cur_lm, scal_l,...
                                                         'B', args.B, 'L', args.L, ...
                                                         'Spin', args.Spin, ...
                                                         'J_min', args.J_min, ...
                                                         'Upsample', args.Upsample,...
                                                         'Reality', args.Reality,...
                                                         'SpinLowered', args.SpinLowered, ...
                                                         'SpinLoweredFrom', args.SpinLoweredFrom,...
                                                         'Sampling', args.Sampling );
% Clear array
cur_lm = 0; 
scal_l = 0; 
f_cur_lmn_syn =0; 
f_scal_lm_syn =0;
f_cur_lmn_syn_rotated =0;                     
end