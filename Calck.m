function k = Calck(h, S, m)
    k = ones(size(h));
    for i = 1:length(h)
        if h(i) < 0
            k(i) = sqrt(S(i)) * (1 - (1 - S(i)^(1/m(i)))^(m(i)))^2;
        end
    end
end