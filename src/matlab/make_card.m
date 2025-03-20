function card=make_card(keyword,value)
%MAKE_CARD turns a set of strings into a valid FITS card
%Make keyword field 8 bytes long
lk=length(keyword);
if (lk > 8) & (nargin>1)
	error('Keyword must be less than or equal to 8 characters!')
elseif (lk < 8 )
	keyword=[keyword,setstr(ones(1,8-lk)*32)];
end;

%Deal with both straight keyword and keyword/value pair
if (nargin==1)
	%Keyword without a value
	card=keyword;	
else
	%Key/value pair has an equal sign and space at bytes 9 and 10
	card=[keyword,'= '];

	%Now output the value. The FITS standard wants things to start 
	%in different columns depending on what type of data the
	%value holds, according to the following rules:
	%
	%  Logical: T or F in column 30
	%
	%  Character string: A beginning quote in column 11 and an
	%  ending quote between columns 20 and 80.
	%
	%  Real part of an integer or floating point number: right 
	%  justified, ending in column 30.
	%
	%  Imaginary part: right justified, ending in
	%  column 50, starting after column 30 (NB. I won't bother 
	%  to support an imaginary part in this M-file, and will 
	%  let some radio astronomer who needs it add it if they want).

	if isstr(value)
  	    %Test for logical. If logical it goes in column 30 
		if (length(value)==1) & (strmatch(upper(value),'T') | strmatch(upper(value),'F'))
 			card=[card,setstr(ones(1,19)*32),value];	
		else	
			%Value must be a character string. Pad if less than 8
			%characters long.
			lv=length(value);
		    if (lv > 70)
		   error('Value must be less than 70 characters long!')
		    elseif (lv < 10 )
 	 	   value=[value,setstr(ones(1,8-lv)*32)];
		    end;
			card=[card,'''',value,''''];
		end;	
	else
		%Value must be a number. Convert to a string. Maximum
		%precision is set to 10 digits
		value=num2str(value,10);
		lv=length(value);
	
		%Write this out so it will end on column 30
		card=[card,setstr(ones(1,20-lv)*32),value];	
	end;
end;

%Now pad the output to make it exactly 80 bytes long
card=[card,setstr(ones(1,80-length(card))*32)];
