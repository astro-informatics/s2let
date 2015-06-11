function error_on_tiling = s2let_check_tiling(psi, phi, L, spin, J)

% s2let_check_axisym_tiling - Checks exactness of the tiling.
% -- Spin directional wavelets on the sphere.
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

identity = zeros(1,L);
for l=abs(spin):L-1
	identity(1,l+1) = identity(1,l+1) + 4*pi/(2*l+1) * phi(l+1) * conj(phi(l+1));
end

for j=0:J
	ind = spin*spin + 1;
	for l=abs(spin):L-1
		for m=-l:l
		    identity(1,l+1) = identity(1,l+1) + 8*pi^2/(2*l+1) * psi(ind, j+1) * conj(psi(ind, j+1));
			ind = ind + 1;
		end
	end
end

error_on_tiling = 0;
for l=abs(spin):L-1
    error_on_tiling = error_on_tiling + identity(1,l+1) - 1.0;
end

end
