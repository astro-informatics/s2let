function vec = s2let_mw_arr2vec(arr)

sz = size(arr);

if sz(1) == 2*sz(2)-1  % array 2L-1 x L
    
    L = sz(2);
    vec = zeros(L*(2*L-1),1);
    for t = 1:L
        for p = 1:2*L-1
            vec((t-1)*(2*L-1)+p,1) = arr(p,t);
        end
    end
    
end

if sz(2) == 2*sz(1)-1  % array L x 2L-1
    
    L = sz(1);
    vec = zeros(L*(2*L-1),1);
    for t = 1:L
        for p = 1:2*L-1
            vec((t-1)*(2*L-1)+p,1) = arr(t,p);
        end
    end
    
end

end