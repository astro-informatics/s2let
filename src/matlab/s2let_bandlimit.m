function bl = s2let_bandlimit(j, J_min, B, L)

% s2let_bandlimit
% Return the band-limit of a specific wavelet scale j
%
% Default usage:
%		bl = s2let_bandlimit(j, J_min, B, L)
%	
%   j the scale of interest,
%   J_min the minimal wavelet scale,
%	B is the wavelet parameter,
%	L is the band-limit for the transform.
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

bl = s2let_bandlimit_mex(j, J_min, B, L);

end
