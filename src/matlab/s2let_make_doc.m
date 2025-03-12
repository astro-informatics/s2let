function s2let_make_doc(s2letpath)

% s2let_make_doc
% Generate Matlab documentation
%
% Default usage :
%
%   s2let_make_doc(s2letpath)
%
% s2letpath is the path for the S2LET package (root)
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

cd(s2letpath)
m2html('mfiles', 'src/main/matlab', 'htmldir', 'doc/matlab');

end