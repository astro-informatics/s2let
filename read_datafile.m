% read datafile "f_cur_lmn.dat" containing "n, el, m, f_cur_lmn" 


%formatSpec= ' %d, %d, %d, %d, %6.5e+i%6.5e\n' ;
%open file
%fid= fopen('f_cur_lmn.dat');
%rawData=fscanf(fid, formatSpec);
%fclose(fid);

fid= fopen('f_cur_lmnONLYj0.dat');
%fid= fopen('f_cur_lmnONLYj1.dat');
%fid= fopen('f_cur_lmnONLYj2.dat');
%fid= fopen('f_cur_lmnONLYj3.dat');
%fid= fopen('f_cur_lmnONLYj4.dat');
rawData2=fscanf(fid, '%f, %f',[2 inf])
%fclose(fid2);
complexData=complex(rawData2(1,:),rawData2(2,:))
