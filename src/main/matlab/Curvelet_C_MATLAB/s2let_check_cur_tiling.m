function error_on_cur_tiling = s2let_check_cur_tiling(psi, phi, L, spin, J)

% s2let_check_cur_tiling - Checks exactness of the tiling.
% -- Curvelets on the sphere.
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

% Scaling function
identity = zeros(1,L);
for l=abs(spin):L-1
	identity(1,l+1) = identity(1,l+1) + 4*pi/(2*l+1) * phi(l+1) * conj(phi(l+1));
end


% Curvelet function
for j=0:J
	ind = spin*spin + 1;
	for l=abs(spin):L-1
		for m= -l:l
           identity(1,l+1) = identity(1,l+1) + 8*pi^2/(2*l+1) * psi(ind, j+1) * conj(psi(ind, j+1));
			ind = ind + 1;
		end
	end
end

error_on_cur_tiling = 0;
for l=abs(spin):L-1
    error_on_cur_tiling = error_on_cur_tiling + identity(1,l+1) - 1.0;
end

end
