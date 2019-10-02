function S = CalcS(h, alpha, n, m)
    S = ones(size(h));
    for i = 1:length(h)
        if h(i) < 0
            S(i) = (1 +  (-alpha(i)*h(i))^(n(i))  )^(-m(i));
        end
        assert(~isinf(S(i)),'error=INF')
    end
    
end