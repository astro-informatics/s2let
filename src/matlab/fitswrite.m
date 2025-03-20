function fitswrite(data, filename, parname, parval)
%Author: R. G. Abraham, Institute of Astronomy, Cambridge University
%        abraham@ast.cam.ac.uk
%Improved by Boris Leistedt bo be compatible with MW and HEALPIX maps

[nrow,ncol]=size(data);
header_cards = [make_card('SIMPLE','T');       ...
      make_card('BITPIX',16);       ...
      make_card('NAXIS',0);          ...
      make_card('EXTEND','T');        ...
      make_card('END');];
  
[nrowbis,ncolbis] = size(header_cards);
n_blanks = 36 - rem(nrowbis,36);
blank_line = setstr(ones(1,80)*32);
header_cards2 = [header_cards; repmat(blank_line,n_blanks,1)];

header_cards2 = [header_cards2;       ...
      make_card('XTENSION','BINTABLE');        ...
      make_card('BITPIX',8);        ...
      make_card('NAXIS',2);          ...
      make_card('NAXIS1',4);      ...
      make_card('NAXIS2',ncol);      ...
      make_card('PCOUNT',0);        ...
      make_card('GCOUNT',1);        ...
      make_card('TFIELDS',1);        ...
      make_card('TTYPE1','SIGNAL');        ...
      make_card('TFORM1','1E');        ...
      make_card('TUNIT1',' ');        ...
      make_card('EXTNAME','BINTABLE');        ...
      make_card(parname, parval);        ...
      make_card('END')];

header_record = make_header_record(header_cards2);
%[ncards,dummy]=size(header_cards);
%fprintf(header_record(1,:));

fid=fopen(filename,'W');
fwrite(fid,header_record','char');

% try to figure out if we need to swap bytes. This is
% imperfect as I don't know the endian-ness of each
% architecture, so I'm only looking for ones I know for 
% sure are big-endian.
friend = computer;
if strmatch(friend,'PCWIN')
   bswap = 'b';
elseif strmatch(friend,'LNX86')
   bswap = 'b';   
elseif strmatch(friend,'ALPHA')
   bswap = 'b';
else
   bswap = 'l';
end

fwrite(fid, data, 'single', 0, 'B');
%fwrite(fid,data,'single',bswap);
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hrec=make_header_record(card_matrix)

[nrow,ncol] = size(card_matrix);
n_blanks = 36 - rem(nrow,36);
blank_line = setstr(ones(1,80)*32);
hrec = [card_matrix; repmat(blank_line,n_blanks,1)];

