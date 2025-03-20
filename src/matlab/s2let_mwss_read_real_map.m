function [mwmaparr, L]= s2let_mwss_read_real_map(file)

% s2let_mw_read_real_map
% Read an MW real map from a FITS file
% Default usage :
%
%   [mwmaparr, L]= s2let_mw_read_real_map(file)
%
% file the name of the input FITS file,
% mwmaparr the output signal read from the file,
% L its resolution.
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

datacell = fitsread(file,'binarytable');
data = datacell{1};
sz = size(data);

mwmap = [];
for col = 1:sz(1)
    mwmap = [mwmap data(col,:)];
end

mwmaparr = s2let_mwss_vec2arr(mwmap);
sz = size(mwmaparr);
L = sz(1)-1;

end
