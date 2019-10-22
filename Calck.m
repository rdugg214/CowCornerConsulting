function k = Calck(hvec, Svec, mvec)
    k = ones(size(hvec));
%     inner = S(locs).*(1./m(locs));
%     innera = (1-inner);
%     innerb = (1-innera).^(m(locs));
%   k(locs) = sqrt(S(locs)) .* innerb.^2;
%  k(locs) = sqrt(S(locs)) .* (1 - (1 - S(locs).^(1./m(locs))).^(m(locs))).^2;

for i = 1:length(hvec)
    S = Svec(i);
    h = hvec(i);
    m = mvec(i);
    if h<0
%         k(i) = sqrt(S(i)) * (1-((1-S(i)^(1/m(i))) )^(m(i)))^2;
            k(i) = sqrt(S) *  ((1-((1- (S^(1/m)))^m))^2);
    end

end
end