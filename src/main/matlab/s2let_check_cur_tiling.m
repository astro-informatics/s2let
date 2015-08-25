function error_on_cur_tiling = s2let_check_cur_tiling(cur_lm, scal_l, L, Spin, J, J_min)
% -----------------------------------------------------------
% s2let_check_cur_tiling 
% - Checks whether the admissibity condition is satisfied. 
% -----------------------------------------------------------
% Log:
% -  constructed by Jennifer Y H Chan on 5th June 2015
% -----------------------------------------------------------
% S2LET package to perform wavelets transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details
% -----------------------------------------------------------

% Scaling function
identity = zeros(1,L);
for el=abs(Spin):L-1
    identity(1,el+1) = identity(1,el+1)+4.*pi/(2*el+1)*scal_l(el^2+el+1,1)*conj(scal_l(el^2+el+1,1)); 
end


% Curvelet functions
for j= J_min: J 
	ind = Spin*Spin + 1;
	for el=abs(Spin):L-1
		for m= -el:el 
          identity(1,el+1) = identity(1,el+1)+(8.*pi^2/(2*el+1))* cur_lm{j-J_min+1}(ind)*conj(cur_lm{j-J_min+1}(ind));
		  ind = ind + 1;
		end
	end
end

error_on_cur_tiling = 0;
for el=abs(Spin):L-1
    error_on_cur_tiling = error_on_cur_tiling+(identity(1,el+1) - 1.0);
end


end

