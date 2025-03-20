function arr = s2let_mwss_vec2arr(vec)

len = length(vec(:));

delta = sqrt(1 + 2*len);
L = ( delta - 1 ) / 2;

arr = zeros(L+1, 2*L);

for t = 1:L+1
    for p = 1:2*L
        arr(t,p) = vec((t-1)*2*L+p);
    end
end


end
