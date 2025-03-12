function error_on_axisym_tiling = s2let_check_axisym_tiling(kappa, kappa0, L, J)

% s2let_check_axisym_tiling - Checks exactness of the tiling.
% -- Axisymmetric wavelets on the sphere.
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

identity = kappa0.^2;
for j=0:J
    identity(1,:) = identity(1,:) + kappa(j+1,:).^2;
end

error_on_axisym_tiling = 0;
for l=1:L
    error_on_axisym_tiling = error_on_axisym_tiling + identity(1,l) - 1.0;
end

end
