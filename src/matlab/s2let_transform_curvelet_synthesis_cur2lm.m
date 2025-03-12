function flm_rec = s2let_transform_curvelet_synthesis_cur2lm(f_cur, f_scal,  varargin)

% s2let_transform_curvelet_synthesis_cur2lm
% Compute curvelet transform:  
% input in curvelet space (i.e. hamonic to Wigner space via SO3_forward_direct)
% output in harmonic space (i.e. Wigner to harmonic space via synthesis_lmn2lm) .
%
% Default usage :
%
% flm_rec = s2let_transform_curvelet_synthesis_lm2cur(f_cur, f_scal,  <options>)
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
%  'Reality'         = { false      [do not assume corresponding signal f real (default)],
%                        true       [assume f real (improves performance)] }
%  'Upsample'        = { false      [multiresolution algorithm (default)],
%                        true       [full resolution curvelets] }
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

len = size(f_cur);
temp = f_cur{len};
sz = size(temp);
Lguessed = sz(2);

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
% Signal synthesis: (Transform to lmn space, then reconstruct the signal in harmonic space)
% ----------------- 
% Scaling function contribution:
if (args.Upsample == 0)  
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
    % Nj = orientational band-limit at j-th scale
    Nj = band_limit;
    if (args.Upsample == 0)  
        f_cur_lmn_syn{j-args.J_min+1} = so3_forward_direct(f_cur{j-args.J_min+1} , band_limit, Nj, ...
                                                    'Reality', args.Reality, 'Sampling', args.Sampling);
    else  
        f_cur_lmn_syn{j-args.J_min+1} = so3_forward_direct(f_cur{j-args.J_min+1} , args.L, Nj, ...
                                                    'Reality', args.Reality, 'Sampling', args.Sampling);
    end
end

% Rotate the Wigner coefficients f_cur_lmn (such the curvelets centered at the North pole)
% Exploit the property of curvelets that cur_ln = (cur_ll)*(delta_ln) 
% such that cur_lmk_rotated = cur_lml*conj(Dlkl(evaluated at the desired rotation angle))
% ---------------
% Define Euler angles
% (for rotating the curvelets to the north pole):
% ---------------
alpha = pi ;
gamma = 0 ;

for j = args.J_min:J,
    band_limit = min([ s2let_bandlimit(j,args.J_min,args.B,args.L) args.L ]);
    Nj = band_limit;
% ---------------
% Precompute Wigner small-d functions
% denoted here as d (in the paper: d_lmn for all el, m, n evaluated at beta).
% They are indexed d(el,m,n).
% Alpha and gamma are the other two rotation angles.
% ---------------
    if (args.Upsample ~= 0)
        beta = pi-acos(-args.Spin/args.B^j);
        d = zeros(args.L, 2*args.L-1, 2*args.L-1);
        d(1,:,:) = ssht_dl(squeeze(d(1,:,:)), args.L, 0, beta);
        for el = 1:args.L-1
            d(el+1,:,:) = ssht_dl(squeeze(d(el,:,:)), args.L, el, beta);
        end
    else
        beta = pi-acos(-args.Spin/args.B^j);
        d = zeros(band_limit, 2*band_limit-1, 2*band_limit-1);
        d(1,:,:) = ssht_dl(squeeze(d(1,:,:)), band_limit, 0, beta);
        for el = 1:band_limit-1
           d(el+1,:,:) = ssht_dl(squeeze(d(el,:,:)), band_limit, el, beta);
        end
    end
  % for the case SO3_STORAGE_PADDED:
    if (args.Reality == 0)  
        if (args.Upsample ~= 0) 
            f_cur_lmn_syn_rotated{j-args.J_min+1} = zeros((2*Nj-1)*args.L^2,1); 
        else 
            f_cur_lmn_syn_rotated{j-args.J_min+1} = zeros((2*Nj-1)*band_limit^2,1);
        end 
        for el = abs(args.Spin):(band_limit-1) 
            for m = -el:el
                if (args.Upsample == 0)  
                         ind_lml = so3_elmn2ind(el,m,el,band_limit,Nj);
                         ind_lmnl = so3_elmn2ind(el,m,-el,band_limit,Nj);
                else
                         ind_lml = so3_elmn2ind(el,m,el,args.L,Nj);
                         ind_lmnl = so3_elmn2ind(el,m,-el,args.L,Nj);
                end 
                en_max = min(el, Nj-1); 
                for en = -en_max:en_max
                     %  Dlmn = exp(-1i*m*alpha) * d(el+1,m+L,n+L) * exp(-1i*n*gamma);
                    if (args.Upsample ~= 0)
                        Dlln = exp(-1i*el*alpha) * d(el+1,el+args.L,en+args.L) * exp(-1i*en*gamma);
                        Dl_nl_n = exp(-1i*-el*alpha) * d(el+1,-el+args.L,en+args.L) * exp(-1i*en*gamma);
                        ind_lmn = so3_elmn2ind(el,m,en,args.L,Nj);
                    else
                        Dlln = exp(-1i*el*alpha) * d(el+1,el+band_limit,en+band_limit) * exp(-1i*en*gamma);
                        Dl_nl_n = exp(-1i*-el*alpha) * d(el+1,-el+band_limit,en+band_limit) * exp(-1i*en*gamma);
                        ind_lmn = so3_elmn2ind(el,m,en,band_limit,Nj);
                     end  
                     f_cur_lmn_syn_rotated{j-args.J_min+1}(ind_lml)=  f_cur_lmn_syn_rotated{j-args.J_min+1}(ind_lml)+ ...
                                                                      conj(Dlln) * f_cur_lmn_syn{j-args.J_min+1}(ind_lmn); 
                     f_cur_lmn_syn_rotated{j-args.J_min+1}(ind_lmnl)= f_cur_lmn_syn_rotated{j-args.J_min+1}(ind_lmnl) + ...
                                                                      conj(Dl_nl_n) * f_cur_lmn_syn{j-args.J_min+1}(ind_lmn); 
                end % end n-loop 
            end  % end m-loop 
        end % end el-loop 
    else % real setting: 
        if (args.Upsample ~= 0) 
            f_cur_lmn_syn_rotated{j-args.J_min+1} = zeros(Nj*args.L^2,1);  
        else
            f_cur_lmn_syn_rotated{j-args.J_min+1} = zeros(Nj*band_limit^2,1);  
        end
        for el = 0:(band_limit-1) 
            for m = -el:el
                % (n=0) terms
                if (args.Upsample == 0)  
                    ind_lml = so3_elmn2ind(el,m,el,band_limit,Nj, 'Reality', args.Reality);
                    ind_lmnzero = so3_elmn2ind(el,m,0,band_limit,Nj, 'Reality', args.Reality);
                    Dl_l_nzero = exp(-1i*el*alpha) * d(el+1,el+band_limit,0+band_limit) * exp(-1i*0*gamma);
                else
                    ind_lml = so3_elmn2ind(el,m,el,args.L,Nj, 'Reality', args.Reality);
                    ind_lmnzero = so3_elmn2ind(el,m,0,args.L,Nj, 'Reality', args.Reality);
                    Dl_l_nzero = exp(-1i*el*alpha) * d(el+1,el+args.L,0+args.L) * exp(-1i*0*gamma);
                end
                f_cur_lmn_syn_rotated{j-args.J_min+1}(ind_lml)=  conj(Dl_l_nzero)*f_cur_lmn_syn{j-args.J_min+1}(ind_lmnzero);                           
                % (n> 0) terms
                en_max = min(el, Nj-1); 
                for en = 1:en_max
                    if (args.Upsample ~= 0)
                        Dl_l_n = exp(-1i*el*alpha) * d(el+1,el+args.L,en+args.L) * exp(-1i*en*gamma);
                        Dl_l_nn = exp(-1i*el*alpha) * d(el+1,el+args.L,-en+args.L) * exp(-1i*-en*gamma);
                        ind_lmn = so3_elmn2ind(el,m,en,args.L,Nj, 'Reality', args.Reality);
                        ind_l_nm_n = so3_elmn2ind(el,-m,en,args.L,Nj, 'Reality', args.Reality);
                    else
                        Dl_l_n = exp(-1i*el*alpha) * d(el+1,el+band_limit,en+band_limit) * exp(-1i*en*gamma);
                        Dl_l_nn = exp(-1i*el*alpha) * d(el+1,el+band_limit,-en+band_limit) * exp(-1i*-en*gamma);
                        ind_lmn = so3_elmn2ind(el,m,en,band_limit,Nj, 'Reality', args.Reality);
                        ind_l_nm_n = so3_elmn2ind(el,-m,en,band_limit,Nj, 'Reality', args.Reality);
                    end
                    if (mod((m+en),2) == 1)
                           sign = -1; 
                    else 
                           sign = 1; 
                    end 
                    f_cur_lmn_syn_rotated{j-args.J_min+1}(ind_lml)=  f_cur_lmn_syn_rotated{j-args.J_min+1}(ind_lml)+ ...
                                                                     conj(Dl_l_n) * f_cur_lmn_syn{j-args.J_min+1}(ind_lmn)+ ... 
                                                                     sign*conj(Dl_l_nn)*conj(f_cur_lmn_syn{j-args.J_min+1}(ind_l_nm_n)); 
                end % end en-loop 
            end  % end m-loop 
        end % end el-loop          
    end % end if Reality loop  
end %end j-loop 

% Reconstruct the signals in harmonic space and perform Wigner transform 
flm_rec = s2let_transform_curvelet_synthesis_lmn2lm(f_cur_lmn_syn_rotated, f_scal_lm_syn, cur_lm, scal_l,...
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