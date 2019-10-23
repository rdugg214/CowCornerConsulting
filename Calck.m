function k = Calck(h, S, m,x,z,dx,dz,hetgen)
    k = ones(size(h));
%     inner = Svec(i)(locs).*(1./mvec(i)(locs));
%     innera = (1-inner);
%     innerb = (1-innera).^(mvec(i)(locs));
%   k(locs) = sqrt(Svec(i)(locs)) .* innerb.^2;
%  k(locs) = sqrt(Svec(i)(locs)) .* (1 - (1 - Svec(i)(locs).^(1./mvec(i)(locs))).^(mvec(i)(locs))).^2;

for i = 1:length(h)
    if h(i)<0 && hetgen.boundary(i) && (z(i) == 100 || z(i) == 0 || x(i) == 0 || x(i) == 500)
    k(i) = sqrt(S(i)) *  ((1-((1- (S(i)^(1/m(i))))^m(i)))^2);
    elseif h(i)<0 && hetgen.boundary(i)
        for j = 1:4
    k(i) = k(i) + (sqrt(S(i)) *  ((1-((1- (S(i)^(1/m(i+hetgen.xcos(j) + ...
        hetgen.zcos(j)))))^m(i+hetgen.xcos(j) +  hetgen.zcos(j))))^2))...
        /((1/4) *dx(i+hetgen.xnos(j))*dz(i+hetgen.znos(j)));
        end
    elseif h(i)<0
%         k(i) = sqrt(Svec(i)(i)) * (1-((1-Svec(i)(i)^(1/mvec(i)(i))) )^(mvec(i)(i)))^2;
            k(i) = sqrt(S(i)) *  ((1-((1- (S(i)^(1/m(i))))^m(i)))^2);
    end

end
end