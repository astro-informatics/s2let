function vec = s2let_mwss_arr2vec(arr)

sz = size(arr);

if sz(1) == 2*(sz(2)-1)  % array 2L x L+1

    L = sz(2)-1;
    vec = zeros((L+1)*2*L,1);
    for t = 1:L+1
        for p = 1:2*L
            vec((t-1)*2*L+p,1) = arr(p,t);
        end
    end

end

if sz(2) == 2*(sz(1)-1)  % array L+1 x 2L

    L = sz(1)-1;
    vec = zeros((L+1)*2*L,1);
    for t = 1:L+1
        for p = 1:2*L
            vec((t-1)*2*L+p,1) = arr(t,p);
        end
    end

end

end
