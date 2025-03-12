% s2let_hpxtest
% Run test for HEALPIX data.
% The wavelet transforms will not be exact but should
% at least reconstruct the input maps with correct precision.
%
% S2LET package to perform Wavelets on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

clear all;
close all;

setenv('HEALPIX','/Users/bl/software/Healpix_2.20a')
setenv('HEALPIXDATA','/Users/bl/software/Healpix_2.20a/data')

% Main parameters
L = 4;
nside = 2;
B = 2;
J_min = 1;
J = s2let_jmax(L, B);
npix = 12*nside*nside;
Spin = 2;


TT_l = rand(L,1);
EE_l = rand(L,1);
TE_l = rand(L,1);
BB_l = rand(L,1);
for l = 0:Spin-1
	TT_l(l+1) = 0;
	EE_l(l+1) = 0;
	TE_l(l+1) = 0;
	BB_l(l+1) = 0;
end
[T_lm, E_lm, B_lm] = ebsep_sim_maps(TT_l, EE_l, TE_l, BB_l, L);
[QpU_lm] = ebsep_compute_QpU_fullsky(E_lm, B_lm);
[QmU_lm] = ebsep_compute_QmU_fullsky(E_lm, B_lm);


disp('Perform spin spherical harmonic decomposition with default parameters')
[fQ, fU] = s2let_hpx_alm2map_spin(E_lm, B_lm, nside, 'L', L);
file = 'mapsQUtemp.fits';
s2let_hpx_write_real_spin_maps(fQ, fU, file);
[fQ_bis, fU_bis, nside]= s2let_hpx_read_real_spin_maps(file);
stop
[E_lm_rec, B_lm_rec] = s2let_hpx_map2alm_spin(fQ_bis, fU_bis, 'L', L);
default = max(max(abs(E_lm-E_lm_rec)))
default = max(max(abs(B_lm-B_lm_rec)))


disp('Perform spin spherical harmonic decomposition with custom parameters')
[fQ, fU] = s2let_hpx_alm2map_spin(E_lm, B_lm, nside, 'L', L, 'Spin', Spin);
[E_lm_rec, B_lm_rec] = s2let_hpx_map2alm_spin(fQ, fU, 'L', L, 'Spin', Spin);
default = max(max(abs(E_lm-E_lm_rec)))
default = max(max(abs(B_lm-B_lm_rec)))




disp('Generates band-limited function')
flm = zeros(L^2,1); 
flm = rand(size(flm)) + sqrt(-1)*rand(size(flm));
flm = 2.*(flm - (1+sqrt(-1))./2);

disp('Constraint on flms to generate real signal')
for el = 0:L-1
   ind = el*el + el + 1;
   flm(ind,1) = real(flm(ind,1));
   for m = 1:el
      ind_pm = el*el + el + m + 1;
      ind_nm = el*el + el - m + 1;
      flm(ind_nm,1) = (-1)^m * conj(flm(ind_pm,1));
   end  
end


disp('Perform spherical harmonic decomposition with default parameters')
f = s2let_hpx_alm2map(flm, nside, 'L', L);
flm_rec = s2let_hpx_map2alm(f, 'L', L);
default = max(max(abs(flm-flm_rec)))


disp('Perform spherical harmonic decomposition with custom parameters')
f = s2let_hpx_alm2map(flm, nside, 'L', L);
flm_rec = s2let_hpx_map2alm(f, 'nside', nside, 'L', L);
custom = max(max(abs(flm-flm_rec)))
 

disp('Perform axisym wavelet transform with default parameters')
[f_wav, f_scal] = s2let_transform_axisym_analysis_hpx(f);
f_rec = s2let_transform_axisym_synthesis_hpx(f_wav, f_scal);
default = max(max(abs(f-f_rec)))


disp('Perform axisym wavelet transform with custom parameters')
[f_wav, f_scal] = s2let_transform_axisym_analysis_hpx(f, 'B', B, 'L', L, 'J_min', J_min);
f_rec = s2let_transform_axisym_synthesis_hpx(f_wav, f_scal, 'B', B, 'L', L, 'J_min', J_min);
custom = max(max(abs(f-f_rec)))
