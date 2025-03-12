function s2let_mwss_write_real_map(f, file)

% s2let_mw_write_real_map
% Write an MW real map to a FITS file
% Default usage :
%
%   s2let_mw_write_real_map(f, file)
%
% f the MW map to be written,
% file the name of the output FITS file.
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

sz = size(f);
L = min([sz(1), sz(2)])-1;

mwmap = [];
for l = 1:L+1
    mwmap = [mwmap f(l,:)];
end

fitswrite(mwmap, file, 'L', L);

end
