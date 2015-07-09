function [f_cur, f_scal] = s2let_transform_spin_curvelet_analysis_lm2cur(flm_init, varargin)

% s2let_transform_analysis_lm2cur
% Compute curvelet transform:
% input in harmonic space (i.e. harmonic to Wigner space via analysis_lm2lmn),
% output in curvelet space (i.e. Wigner space to curvelet space via SO3_inverse).
%
% Default usage :
%
%   [f_cur, f_scal] = s2let_transform_curvelet_analysis_lm2cur(flm_init, <options>)
%
% flm_init is the input field in harmonic space,
% f_cur contains the output curvelet contributions,
% f_scal contains the output scaling contributions,
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

sz = length(flm_init(:));
Lguessed = sqrt(sz);

p = inputParser;
p.addRequired('flm_init', @isnumeric);
p.addParamValue('B', 2, @isnumeric);
p.addParamValue('L', Lguessed, @isnumeric);
p.addParamValue('J_min', 0, @isnumeric);
p.addParamValue('Spin', 0, @isnumeric);
p.addParamValue('Reality', false, @islogical);
p.addParamValue('Upsample', false, @islogical);
p.addParamValue('SpinLowered', false, @islogical);
p.addParamValue('SpinLoweredFrom', 0, @isnumeric);
p.addParamValue('Sampling', 'MW', @ischar);
p.parse(flm_init, varargin{:});
args = p.Results;

% For curvelets, azimuthal/directional band-limit N always equals to L           
N = args.L ;
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
% Signal analysis:
% -----------------
% Decompose the signals using curvelets and the scaling functions
% And perform Wigner transform (lm2lmn -  Call matlab function s2let_transform_analysis_lm2lmn)
[f_cur_lmn, f_scal_lm] = s2let_transform_spin_curvelet_analysis_lm2lmn(flm_init, cur_lm, scal_l,...
                                                                      'B',args.B, 'L', args.L, 'J_min', args.J_min, ...
                                                                      'Spin', args.Spin,'Reality', args.Reality,...
                                                                      'Upsample', args.Upsample, ...
                                                                      'SpinLowered', args.SpinLowered, ...
                                                                      'SpinLoweredFrom',  args.SpinLoweredFrom, ...
                                                                      'Sampling', args.Sampling);
% disp('--check the size of f_cur_lmn--')
% For so3_storage_padded: 
 % (2N-1)*L^2 for complex signals 
 %  N*L^2 for real signals  
%whos f_cur_lmn
%len= length(f_cur_lmn)
%temp = f_cur_lmn{len};
%sz = size(temp)
%disp('--------')


% Curvelets:
% Rotate the Wigner coefficients f_cur_lmn (such the curvelets centred at the North pole)
% Exploit the property of curvelets that cur_ln = (cur_ll)*(delta_ln)
% such that cur_lmk_rotated = cur_lml*conj(Dlkl(evaluated at the desired rotation angle))
% ---------------
% Define Euler angles
% (for rotating the curvelets to the north pole):
% ---------------
alpha = 0; % pi ;
beta = pi/2 ;
gamma = 0 ;
% ---------------
% Precompute Wigner small-d functions
% denoted here as d (in the paper: d_lmn for all el, m, n evaluated at beta).
% They are indexed d(el,m,n).
% Alpha and gamma are the other two rotation angles.
% ---------------
d = zeros(args.L, 2*args.L-1, 2*args.L-1);  % check BAND_LIMIT
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
          f_cur_lmn_rotated{j-args.J_min+1} = zeros((2*Nj-1)*args.L^2,1); 
        else 
          f_cur_lmn_rotated{j-args.J_min+1} = zeros((2*Nj-1)*band_limit^2,1);
        end 
        % ind_lmk = 0;
        % ind_lml = 0;
        % ind_lmnl= 0;
        for el = abs(args.Spin):(band_limit-1) %0: args.L-1  % max(abs(args.Spin),abs(en))
            for m = -el:el
                if (args.Upsample == 0)  %false => multi-resolution
                   ind_lml = so3_elmn2ind(el,m,el,band_limit,Nj);
                   ind_lmnl = so3_elmn2ind(el,m,-el,band_limit,Nj);
                else
                   ind_lml = so3_elmn2ind(el,m,el,args.L,Nj);
                   ind_lmnl = so3_elmn2ind(el,m,-el,args.L,Nj);
                end % end the if-loop for upsample
                en_max = min(el, Nj-1); 
                for k = -en_max:en_max 
                        % Dlmn = exp(-1i*m*alpha) * d(el+1,m+L,n+L) * exp(-1i*n*gamma);
                        Dlkl = exp(-1i*k*alpha) * d(el+1,k+args.L,el+args.L) * exp(-1i*el*gamma);  % check BAND_LIMIT
                        Dlknl = exp(-1i*k*alpha) * d(el+1,k+args.L,-el+args.L) * exp(-1i*(-el)*gamma);
                        if (args.Upsample == 0)  %false => multi-resolution
                           ind_lmk = so3_elmn2ind(el,m,k,band_limit,Nj);
                        else
                           ind_lmk = so3_elmn2ind(el,m,k,args.L,Nj);
                        end % end the if-loop for upsample
                        f_cur_lmn_rotated{j-args.J_min+1}(ind_lmk)= conj(Dlkl) * f_cur_lmn{j-args.J_min+1}(ind_lml)+ ...
                                                                    conj(Dlknl)* f_cur_lmn{j-args.J_min+1}(ind_lmnl);                                        
                end % end k-loop
            end % end m-loop 
        end % end el-loop   
    else %i.e. real signals
        if (args.Upsample ~= 0) 
          f_cur_lmn_rotated{j-args.J_min+1} = zeros(Nj*args.L^2,1);  
        else
          f_cur_lmn_rotated{j-args.J_min+1} = zeros(Nj*band_limit^2,1);  
        end
        for el = 0:(band_limit-1) 
            for m = -el:el
                if (args.Upsample == 0)  %false => multi-resolution
                   ind_lml = so3_elmn2ind(el,m,el,band_limit,Nj);
                else
                   ind_lml = so3_elmn2ind(el,m,el,args.L,Nj);
                end % end the if-loop for upsample
                % (k=0) terms
                Dlkzerol = exp(-1i*0*alpha) * d(el+1,0+args.L,el+args.L) * exp(-1i*el*gamma);
                if (args.Upsample == 0)  %false => multi-resolution
                         ind_lmkzero = so3_elmn2ind(el,m,0,band_limit,Nj, 'Reality', args.Reality);
                else
                         ind_lmkzero = so3_elmn2ind(el,m,0,args.L,Nj,'Reality', args.Reality);
                end % end the if-loop for upsample   
                f_cur_lmn_rotated{j-args.J_min+1}(ind_lmkzero)= conj(Dlkzerol) * f_cur_lmn{j-args.J_min+1}(ind_lml); 
                % (k~=0) terms
                en_max = min(el, Nj-1); 
                for k = 1:en_max  %1-mod(en_max,2): en_max
                     Dlkl = exp(-1i*k*alpha) * d(el+1,k+args.L,el+args.L) * exp(-1i*el*gamma);
                     if (args.Upsample == 0)  %false => multi-resolution
                         ind_lmk = so3_elmn2ind(el,m,k,band_limit,Nj, 'Reality', args.Reality);
                     else
                         ind_lmk = so3_elmn2ind(el,m,k,args.L,Nj,'Reality', args.Reality);
                     end % end the if-loop for upsample
                         f_cur_lmn_rotated{j-args.J_min+1}(ind_lmk)= conj(Dlkl) * f_cur_lmn{j-args.J_min+1}(ind_lml); 
                end  % end k-loop
            end % end m-loop
        end % end el-loop
    end % end if (reality)-loop
end %end j-loop

%Debugging
%{
% Curvelets:
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
  ind_lml = 0;
  ind_lmnl =0; 
  ind_lmn = 0;
  for el = abs(args.Spin):(band_limit-1) %0: args.L-1  % max(abs(args.Spin),abs(en))
      for m = -el:el
          en_max = min(el, Nj-1); 
          if (args.Reality == 0) %i.e. false (default) => complex
                 if (args.Upsample == 0)  %false => multi-resolution
                         ind_lml = so3_elmn2ind(el,m,el,band_limit,Nj);
                         ind_lmnl = so3_elmn2ind(el,m,-el,band_limit,Nj);
                 else
                         ind_lml = so3_elmn2ind(el,m,el,args.L,Nj);
                         ind_lmnl = so3_elmn2ind(el,m,-el,args.L,Nj);
                 end % end the if-loop for upsample
                 for en = -en_max:en_max
                     %  Dlmn = exp(-1i*m*alpha) * d(el+1,m+L,n+L) * exp(-1i*n*gamma);
                      Dlln = exp(-1i*el*alpha) * d(el+1,el+args.L,en+args.L) * exp(-1i*en*gamma);
                      Dl_nl_n = exp(-1i*-el*alpha) * d(el+1,-el+args.L,en+args.L) * exp(-1i*en*gamma);
                     if (args.Upsample == 0)  %false => multi-resolution
                         ind_lmn = so3_elmn2ind(el,m,en,band_limit,Nj);
                     else
                         ind_lmn = so3_elmn2ind(el,m,en,args.L,Nj);
                     end % end the if-loop for upsample
                     f_cur_lmn_syn_rotated{j-args.J_min+1}(ind_lml)=  f_cur_lmn_syn_rotated{j-args.J_min+1}(ind_lml) + ...
                                                                      conj(Dlln) * f_cur_lmn_rotated{j-args.J_min+1}(ind_lmn); 
                     f_cur_lmn_syn_rotated{j-args.J_min+1}(ind_lmnl)=  f_cur_lmn_syn_rotated{j-args.J_min+1}(ind_lmnl) + ...
                                                                       conj(Dl_nl_n) * f_cur_lmn_rotated{j-args.J_min+1}(ind_lmn); 
                 end % end n-loop 
          end % end if-loop for reality option 
      end % end m-loop 
  end % end el-loop                    
end
%}


%whos f_cur_lmn_rotated
%len= length(f_cur_lmn_rotated)
%temp = f_cur_lmn{len};
%sz = size(temp)

% -----------------                                                     
% Transform to pixel space:
% -----------------
% Scaling functions:
if (args.Upsample == 0)  %false => multi-resolution 
     band_limit = min([s2let_bandlimit(args.J_min-1,args.J_min,args.B,args.L) args.L ]);
else
     band_limit = args.L ;
end
f_scal = ssht_inverse(f_scal_lm, band_limit,  ...
                      'Method', args.Sampling, ...
                      'Spin', 0, ...
                      'Reality', args.Reality);          
                  



% Rotated-curvelets contributions:
% Then compute the real space curvelet coefficients
for j = args.J_min:J,
    band_limit = min([s2let_bandlimit(j,args.J_min,args.B,args.L) args.L ]);
    Nj = band_limit; 
    if (args.Upsample == 0)  % i.e. false (default) => multi-resolution
        f_cur{j-args.J_min+1} = so3_inverse(f_cur_lmn_rotated{j-args.J_min+1}, band_limit, Nj, ...
                                            'Sampling', args.Sampling, 'Reality', args.Reality) ;
    else
        f_cur{j-args.J_min+1} = so3_inverse(f_cur_lmn_rotated{j-args.J_min+1}, args.L, Nj, ...
                                            'Sampling', args.Sampling, 'Reality', args.Reality) ;
       %For debugging:                                 
       % size(f_cur{j-args.J_min+1}) 
    end
end
% For debugging: 
%disp('--check the size of f_cur--')
%whos f_cur
%len= length(f_cur)
%temp = f_cur{len};
%sz = size(temp)
%disp('--------')
% For L=N=16, real data: sz = 31  16 31  
%{
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
%}


% Clear array
cur_lm = 0; 
scal_l = 0; 
f_cur_lmn =0; 
f_scal_lm =0;

end