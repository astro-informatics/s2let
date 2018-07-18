function [f_cur, f_scal] = s2let_transform_curvelet_analysis_lm2cur(flm_init, varargin)

% s2let_transform_curvelet_analysis_lm2cur
% Compute curvelet transform:
% input in harmonic space  (i.e. harmonic to Wigner space via analysis_lm2lmn),
% output in curvelet space (i.e. Wigner space to curvelet space via SO3_inverse_direct).
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
%
% -----------------------------------------------------------
% S2LET package to perform Wavelet Transform on the Sphere.
% Copyright (C) 2012-2016  Boris Leistedt, Jennifer Chan & Jason McEwen
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

J = s2let_jmax(args.L, args.B);

% ---------------
% Tile curvelets:
% ---------------
% Call curvelet- and scaling-function- generating functions
[cur_lm scal_l] = s2let_curvelet_tiling(args.B, args.L, args.J_min,...
                                        'Spin', args.Spin,...
                                        'SpinLowered', args.SpinLowered, ...
                                        'SpinLoweredFrom', args.SpinLoweredFrom);

% -----------------
% Signal analysis:
% -----------------
% Decompose the signals using curvelets and the scaling functions 
% Then perform Wigner transform (calling matlab function s2let_transform_analysis_lm2lmn)
[f_cur_lmn, f_scal_lm] = s2let_transform_curvelet_analysis_lm2lmn(flm_init, cur_lm, scal_l,...
                                                                  'B',args.B, 'L', args.L, 'J_min', args.J_min, ...
                                                                  'Spin', args.Spin,'Reality', args.Reality,...
                                                                  'Upsample', args.Upsample, ...
                                                                  'SpinLowered', args.SpinLowered, ...
                                                                  'SpinLoweredFrom',  args.SpinLoweredFrom, ...
                                                                  'Sampling', args.Sampling);

% Curvelet contribution:
% Rotate the Wigner coefficients f_cur_lmn (such the curvelets centered at the North pole)
% Exploit the property of curvelets that cur_ln = (cur_ll)*(delta_ln) 
% such that cur_lmk_rotated = cur_lml*conj(Dlkl(evaluated at the desired rotation angle))
% ---------------
% Define Euler angles for rotation: 
% ---------------
alpha = 0;
gamma = 0 ;

for j = args.J_min:J,
    band_limit = min([ s2let_bandlimit(j,args.J_min,args.B,args.L) args.L ]);
    % Nj = orientational band-limit at j-th scale:
    Nj = band_limit;
% ---------------
% Specify beta and orecompute Wigner small-d functions, denoted here as d (in the paper: d_lmn for all el, m, n evaluated at beta).
% They are indexed d(el,m,n). Alpha and gamma are the other two rotation angles.
% ---------------
    if (args.Upsample ~= 0)
        beta = acos(-args.Spin/args.B^j);
        d = zeros(args.L, 2*args.L-1, 2*args.L-1);
        d(1,:,:) = ssht_dl(squeeze(d(1,:,:)), args.L, 0, beta);
        for el = 1:args.L-1
            d(el+1,:,:) = ssht_dl(squeeze(d(el,:,:)), args.L, el, beta);
        end
    else
        beta = acos(-args.Spin/args.B^j);
        d = zeros(band_limit, 2*band_limit-1, 2*band_limit-1);
        d(1,:,:) = ssht_dl(squeeze(d(1,:,:)), band_limit, 0, beta);
        for el = 1:band_limit-1
            d(el+1,:,:) = ssht_dl(squeeze(d(el,:,:)), band_limit, el, beta);
        end
    end
    % for the case SO3_STORAGE_PADDED:
    if (args.Reality == 0) 
        if (args.Upsample ~= 0) 
          f_cur_lmn_rotated{j-args.J_min+1} = zeros((2*Nj-1)*args.L^2,1); 
        else 
          f_cur_lmn_rotated{j-args.J_min+1} = zeros((2*Nj-1)*band_limit^2,1);
        end 
        for el = abs(args.Spin):(band_limit-1) 
            for m = -el:el
                if (args.Upsample == 0)   
                   ind_lml = so3_elmn2ind(el,m,el,band_limit,Nj);
                   ind_l_m_nl = so3_elmn2ind(el,m,-el,band_limit,Nj);
                else  
                   ind_lml = so3_elmn2ind(el,m,el,args.L,Nj);
                   ind_l_m_nl = so3_elmn2ind(el,m,-el,args.L,Nj);
                end 
                en_max = min(el, Nj-1); 
                for k = -en_max:en_max 
                        % Dlmn = exp(-1i*m*alpha) * d(el+1,m+L,n+L) * exp(-1i*n*gamma);
                    if (args.Upsample ~= 0)
                        Dlkl = exp(-1i*k*alpha) * d(el+1,k+args.L,el+args.L) * exp(-1i*el*gamma);  
                        Dlknl = exp(-1i*k*alpha) * d(el+1,k+args.L,-el+args.L) * exp(-1i*(-el)*gamma);
                        ind_lmk = so3_elmn2ind(el,m,k,args.L,Nj);
                    else
                        Dlkl = exp(-1i*k*alpha) * d(el+1,k+band_limit,el+band_limit) * exp(-1i*el*gamma);
                        Dlknl = exp(-1i*k*alpha) * d(el+1,k+band_limit,-el+band_limit) * exp(-1i*(-el)*gamma);
                        ind_lmk = so3_elmn2ind(el,m,k,band_limit,Nj);
                    end
                    f_cur_lmn_rotated{j-args.J_min+1}(ind_lmk)= conj(Dlkl) * f_cur_lmn{j-args.J_min+1}(ind_lml)+ ...
                                                                    conj(Dlknl)* f_cur_lmn{j-args.J_min+1}(ind_l_m_nl);                                        
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
                if (args.Upsample == 0)  
                   ind_lml = so3_elmn2ind(el,m,el,band_limit,Nj, 'Reality', args.Reality); 
                   ind_l_nm_l = so3_elmn2ind(el,-m,el,band_limit,Nj, 'Reality', args.Reality);
                else
                   ind_lml = so3_elmn2ind(el,m,el,args.L,Nj, 'Reality', args.Reality) ; 
                   ind_l_nm_l = so3_elmn2ind(el,-m,el,args.L,Nj, 'Reality', args.Reality) ;
                end 
                if (mod((m+el),2) == 1) 
                    sign = -1; 
                else     
                    sign = 1; 
                end 
                en_max = min(el, Nj-1); 
                for k = 0:en_max
                     if (args.Upsample ~= 0)
                         Dl_k_l = exp(-1i*k*alpha) * d(el+1,k+args.L,el+args.L) * exp(-1i*el*gamma);
                         Dl_k_nl = exp(-1i*k*alpha) * d(el+1,k+args.L,-el+args.L) * exp(-1i*-el*gamma);
                         ind_lmk = so3_elmn2ind(el,m,k,args.L,Nj,'Reality', args.Reality);
                     else
                         Dl_k_l = exp(-1i*k*alpha) * d(el+1,k+band_limit,el+band_limit) * exp(-1i*el*gamma);
                         Dl_k_nl = exp(-1i*k*alpha) * d(el+1,k+band_limit,-el+band_limit) * exp(-1i*-el*gamma);
                         ind_lmk = so3_elmn2ind(el,m,k,band_limit,Nj, 'Reality', args.Reality);
                     end 
                     f_cur_lmn_rotated{j-args.J_min+1}(ind_lmk)= conj(Dl_k_l) * f_cur_lmn{j-args.J_min+1}(ind_lml)+ ...
                                                                 sign*conj(Dl_k_nl)* conj(f_cur_lmn{j-args.J_min+1}(ind_l_nm_l));
                end  % end k-loop
            end % end m-loop
        end % end el-loop
    end % end if (reality)-loop
end %end j-loop


% -----------------                                                     
% Transform to pixel space:
% -----------------
% Scaling functions in real space:
if (args.Upsample == 0)  
     band_limit = min([s2let_bandlimit(args.J_min-1,args.J_min,args.B,args.L) args.L ]);
else
     band_limit = args.L ;
end
f_scal = ssht_inverse(f_scal_lm, band_limit,  ...
                      'Method', args.Sampling, ...
                      'Spin', 0, ...
                      'Reality', args.Reality);          
                  
% Rotated-curvelets contributions:
% Compute the curvelet coefficients in real space
for j = args.J_min:J,
    band_limit = min([s2let_bandlimit(j,args.J_min,args.B,args.L) args.L ]);
    Nj = band_limit; 
    if (args.Upsample == 0)  
        f_cur{j-args.J_min+1} = so3_inverse_direct(f_cur_lmn_rotated{j-args.J_min+1}, band_limit, Nj, ...
                                            'Sampling', args.Sampling, 'Reality', args.Reality) ;
    else
        f_cur{j-args.J_min+1} = so3_inverse_direct(f_cur_lmn_rotated{j-args.J_min+1}, args.L, Nj, ...
                                            'Sampling', args.Sampling, 'Reality', args.Reality) ;
    end
end
% size(f_cur_lmn_rotated{J-args.J_min+1})
% size(f_cur{J-args.J_min+1})


% Clear array
cur_lm = 0; 
scal_l = 0; 
f_cur_lmn =0; 
f_scal_lm =0;

end