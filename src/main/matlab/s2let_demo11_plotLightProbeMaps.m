% s2let_demo11_plotLightProbeMaps
% -----------------------------------------------------------
% Plot the 'light probe' maps from FITS file
% -----------------------------------------------------------
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details
% 
% Modified S2LET package to perform Curvelets on the Sphere.
% -----------------------------------------------------------
% L is the angular band-limit.
% -----------------------------------------------------------
% Log: 
% -  constructed by Jennifer Y H Chan on 11th Aug 2015
% -----------------------------------------------------------
% S2LET package to perform wavelet transform on the Sphere.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details
% -----------------------------------------------------------
 
clear; % all;
% close all ;

% ---------------
% Set the reality and upsample flag for signal analysis:
% ---------------
reality = true;    % true for real data;
upsample = true ;  % true for full resolution plot

% ---------------
% Load data set: 
% ---------------
%filename= '../../../data/heal_grace.fits';
%filename= '../../../data/heal_galileo.fits';
%filename= '../../../data/heal_rnl.fits';
%filename= '../../../data/heal_uffizi.fits';
filename= '../../../data/heal_stpeters.fits';

% ----------
% Check what is in the FITS file
% ----------
info =fitsinfo(filename);
disp(info.Contents);
info.BinaryTable
rowend = info.BinaryTable.Rows;

datacell = fitsread(filename,'binarytable');
data = datacell{1};
sz = size(data);


[data_hpxmap, nside] = s2let_hpx_read_real_map(filename);
whos hpxmap




f_mw = s2let_hpx2mw(data_hpxmap,'nside', 256, 'L',512)
whos data_hpxmap
whos f_mw


L = 32; 
flm_gen= data;      % i.e. cell2mat(flm) -> [768*1024] 
isnumeric(flm_gen)
for row = 1: rowend
   flm_data{row} =   flm_gen(row,:); 
   f_gen{row}  = ssht_inverse(flm_data{row} , L,'Reality', reality);
end
%{
% --------
% Reshape the data from cell to array size appropriate for making MW map
% (i.e vector -> array of size of (L, 2*L-1))
% --------
mwmap = [];
for col = 1:sz(1)
    mwmap = [mwmap data(col,:)];
end

L = sqrt(length(mwmap(:)) /rowend)
mwmaparr= zeros(L, 2*L-1);
for t = 1:L
    for p = 1:2*L-1
        mwmaparr(t,p) = mwmap((t-1)*(2*L-1)+p);
    end
end
size(mwmaparr)
figure
ssht_plot_mollweide(mwmaparr, L, 'Mode', 1)  % Mode 1: plot real part
colorbar;
%v = caxis;
%temp = max(abs(v));
%caxis([-temp temp])
% From the plot: the data set is flm instead of f 
%} 



whos f_gen 


