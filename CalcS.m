function S = CalcS(h, alpha, n, m)
    S = ones(size(h));
     locs = h<0;
     S(locs) = (1 +  (-alpha(locs).*h(locs)).^(n(locs))  ).^(-m(locs));
% for i = 1:length(h)
%     if h(i)<0 && ~hetgen.boundary(i)
%     S(i) = (1+(-alpha(i)*h(i))^(n(i)))^(-m(i));
%     
% end
%   S(locs) = (1+(-alpha(locs).*h(locs)).^n(locs)).^(-m(locs));
%     
end