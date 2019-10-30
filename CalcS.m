function S = CalcS(h, alpha, n, m)
    S = zeros(size(h));
%      locs = h<0;
%      S(locs) = (1 +  (-alpha(locs).*h(locs)).^(n(locs))  ).^(-m(locs));
for i = 1:length(h)
    if h(i)<0
    S(i) = (1+(-alpha(i)*h(i))^(n(i)))^(-m(i));
    else
        S(i) = 1;
    end
end
%   S(locs) = (1+(-alpha(locs).*h(locs)).^n(locs)).^(-m(locs));
%     
end