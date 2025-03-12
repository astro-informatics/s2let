function arr = s2let_mw_vec2arr(vec)

len = length(vec(:));

delta = sqrt(1 + 8*len);
L = ( 1 + delta ) / 4;

arr = zeros(L, 2*L-1);

for t = 1:L
    for p = 1:2*L-1
        arr(t,p) = vec((t-1)*(2*L-1)+p);
    end
end


end
