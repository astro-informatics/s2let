function [healpixmap, nside]= s2let_hpx_read_real_map(file)

% s2let_hpx_read_real_map 
% Read a Healpix real map from a FITS file
% Default usage :
%
%   [healpixmap, L]= s2let_hpx_read_real_map(file)
%
% file the name of the input FITS file,
% healpixmap the output signal read from the file,
% nside its resolution.
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

datacell = fitsread(file,'binarytable');
data = datacell{1};
sz = size(data);

healpixmap = [];
for col = 1:sz(1)
    healpixmap = [healpixmap data(col,:)];
end

nside = sqrt((sz(1) * sz(2)) / 12);

end
