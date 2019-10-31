function S = CalcS(h, alpha, n, m,x,z,dx,dz,hetgen)
%% Homogenous    
S = ones(size(h));
     locs = h<0;
     S(locs) = (1 +  (-alpha(locs).*h(locs)).^(n(locs))  ).^(-m(locs));
%% Heterogenous
% for i = 1:length(h)
%     if h(i)<0 && hetgen.boundary(i) && (z(i) == 100 || z(i) == 0 || x(i) == 0 || x(i) == 500)
%      S(i) = (1+(-alpha(i)*h(i))^(n(i)))^(-m(i));
%     elseif h(i)<0 && hetgen.boundary(i)
%     S(i) = 0;
%         for j = 1:4
% %     S(i) = S(i) + (sqrt(S(i)) *  ((1-((1- (S(i)^(1/m(i+hetgen.xcos(j) + ...
% %         hetgen.zcos(j)))))^m(i+hetgen.xcos(j) +  hetgen.zcos(j))))^2))...
% %         /((1/4) *dx(i+hetgen.xnos(j))*dz(i+hetgen.znos(j)));
%         S(i) = S(i) +((1+(-alpha(i+hetgen.xcos(j) + ...
%         hetgen.zcos(j))*h(i))^(n(i+hetgen.xcos(j) + ...
%         hetgen.zcos(j))))^(-m(i+hetgen.xcos(j) + ...
%         hetgen.zcos(j))))/((1/4) *dx(i+hetgen.xnos(j))*dz(i+hetgen.znos(j)));
%         end
%     elseif h(i)<0
%     S(i) = (1+(-alpha(i)*h(i))^(n(i)))^(-m(i));
%     else
%     S(i) = 1;
%     end
% end
end