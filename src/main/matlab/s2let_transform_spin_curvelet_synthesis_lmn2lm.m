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

sz = length(scal_l(:));
Lguessed = sqrt(sz);

p = inputParser;
p.addRequired('f_cur_lmn', @iscell);
p.addRequired('f_scal_lm', @isnumeric);
p.addRequired('cur_lm', @iscell);
p.addRequired('scal_l', @isnumeric);
p.addParamValue('B', 2, @isnumeric);
p.addParamValue('L', Lguessed, @isnumeric);   
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
J = s2let_jmax(args.L, args.B);

% ---------------
% Signal synthesis: 
% ---------------
% Generate flmn from flm of the complex signal
% (c.f. C function: s2let_synthesis_cur_lmn2lm(flm, f_cur_lmn, f_scal_lm, cur_lm, scal_l, parameters);)
% -----------------
% Curvelet contribution:
% -----------------
% disp('syn_lmn2lm: Curvelet synthesis of complex signals from Wigner to harmonic space (flm_cur_syn from flmn_syn)');
flm_rec = zeros(args.L^2,1);
for j = args.J_min:J,
    ind_ln =0;
    ind_lnn =0;
    ind_lm=0;
    ind_lmn=0;
    band_limit = min([ s2let_bandlimit(j,args.J_min,args.B,args.L) args.L ]);
    Nj = band_limit;
    %if (args.Upsample ~= 0)  %true => full-resolution
    %    band_limit = args.L;
    %end
    if (args.Reality == 0) %i.e. false (default) => complex
        for en = -Nj+1:Nj-1,
            for el = max(abs(args.Spin),abs(en)):band_limit-1,
                ind_ln = ssht_elm2ind(el, en);
                psi = (cur_lm{j-args.J_min+1}(ind_ln));
                for m = -el:el,
                    ind_lm = ssht_elm2ind(el, m);
                    if (args.Upsample ~= 0)  %true => full-resolution
                        ind_lmn = so3_elmn2ind(el,m,en,args.L,Nj);
                    else  %false => multi-resolution
                        ind_lmn = so3_elmn2ind(el,m,en,band_limit,Nj);
                    end % end the if-loop for upsample
                    flm_rec(ind_lm)= flm_rec(ind_lm)+ f_cur_lmn{j-args.J_min+1}(ind_lmn)* psi;
                end
            end
        end
    else % i.e.(args.Reality == 1) %i.e. true => real
        for el = 0:band_limit-1,     
            for m = -el:el,
                ind_lm = ssht_elm2ind(el, m);
                % contribution from n=0 terms:
                 ind_lnzero = ssht_elm2ind(el, 0);
                 if (args.Upsample ~= 0) 
                     ind_lmnzero = so3_elmn2ind(el,m,0,args.L,Nj,'Reality', args.Reality);
                 else
                     ind_lmnzero = so3_elmn2ind(el,m,0,band_limit,Nj,'Reality', args.Reality);
                 end
                 psizero= cur_lm{j-args.J_min+1}(ind_lnzero);  
                 flm_rec(ind_lm)=  flm_rec(ind_lm)+ f_cur_lmn{j-args.J_min+1}(ind_lmnzero)* psizero; 
                % for n ~= 0 terms: 
                 if (args.Upsample ~= 0) 
                     ind_lml = so3_elmn2ind(el,m,el,args.L,Nj,'Reality', args.Reality);
                     ind_l_nm_l = so3_elmn2ind(el,-m,el,args.L,Nj,'Reality', args.Reality);  
                 else
                     ind_lml = so3_elmn2ind(el,m,el,band_limit,Nj,'Reality', args.Reality);
                     ind_l_nm_l = so3_elmn2ind(el,-m,el,band_limit,Nj,'Reality', args.Reality);     
                 end
                 if (mod((m+el),2) == 1)    
                     sign = -1; 
                 else      
                     sign = 1; 
                 end 
                 for en = 1:Nj-1,   %-mod(Nj,2)
                     if (el >= en)
                     % contribution from positive n terms 
                      % i.e. Sum_over_lmn((cur_l_n)*f_cur_l_m_n))
                       ind_ln = ssht_elm2ind(el, en); 
                       psi = cur_lm{j-args.J_min+1}(ind_ln);  
                       flm_rec(ind_lm)= flm_rec(ind_lm)+ f_cur_lmn{j-args.J_min+1}(ind_lml)* psi;
                      % contribution from negative n terms: 
                      % i.e. Sum_over_lmn((cur_l_neg-n)*f_cur_l_m_neg-n)) =Sum_over_lmn((cur_l_neg-n)*(-1^(m+n))*f_cur_l_neg-m_n))
                       ind_l_nn = ssht_elm2ind(el, -en); 
                       npsi = cur_lm{j-args.J_min+1}(ind_l_nn);
                       flm_rec(ind_lm)= flm_rec(ind_lm)+ sign*conj(f_cur_lmn{j-args.J_min+1}(ind_l_nm_l))* npsi; 
                     end   % end if (el>=en) loop 
                 end % end en-loop for Reality Option
            end % end em-loop for Reality Option
        end % end el-loop for Reality Option      
    end  % end if-loop for Reality Option
end % end j-loop

% -----------------
% Adding the scaling function contribution: 
% -----------------
% disp('syn_lmn2lm:  Compute flm_scal_syn ');
if (args.Upsample ~= 1)  %false => multi-resolution 
    band_limit = min([ s2let_bandlimit(args.J_min-1, args.J_min, args.B,args.L) args.L ]);
else 
    band_limit = args.L ; 
end
lm_ind=0;
if (args.Reality == 0) %i.e. false (default) => complex
    for el = abs(args.Spin): band_limit-1,
        phi = sqrt(4.*pi/(2.*el+1))*scal_l(el^2+el+1,1);      
        for m = -el:el,
            lm_ind=ssht_elm2ind(el, m);
            flm_rec(lm_ind) =  flm_rec(lm_ind)+  f_scal_lm(lm_ind)* phi;
        end
    end 
else   %  i.e. (args.Reality == 1); true => real
    for el = 0 :band_limit-1,
        phi = sqrt(4.*pi/(2.*el+1))*scal_l(el^2+el+1,1);      
        for m = -el:el,
            lm_ind=ssht_elm2ind(el, m);
            flm_rec(lm_ind) =  flm_rec(lm_ind)+  f_scal_lm(lm_ind)* phi;
        end
    end 
end  % end if-loop for Reality Option

end






