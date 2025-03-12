function s2let_hpx_write_real_map(f, file)

% s2let_hpx_rite_real_map 
% Write a real Healpix map to a FITS file
% Default usage :
%
%   s2let_hpx_write_real_map(f, file)
%
% f the Healpix map to be written,
% file the name of the output FITS file.
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

sz = size(f);
szb = max([sz(1), sz(2)]);
nside = floor(sqrt(szb/12.0));

fitswrite(f, file, 'NSIDE', nside);

end