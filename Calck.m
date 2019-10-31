function k = Calck(h, S, m,x,z,dx,dz,hetgen)
%% Homoogenous
    k = ones(size(h));
    locs = h<0;
    k(locs) = sqrt(S(locs)) .* (1 - (1 - S(locs).^(1./m(locs))).^(m(locs))).^2;
%% Heterogenous
% for i = 1:length(h)
%     if h(i)<0 && hetgen.boundary(i) && (z(i) == 100 || z(i) == 0 || x(i) == 0 || x(i) == 500)
%     k(i) = sqrt(S(i)) *  ((1-((1- (S(i)^(1/m(i))))^m(i)))^2);
%     elseif h(i)<0 && hetgen.boundary(i)
%         k(i) = 0;
%         for j = 1:4
%     k(i) = k(i) + (sqrt(S(i)) *  ((1-((1- (S(i)^(1/m(i+hetgen.xcos(j) + ...
%         hetgen.zcos(j)))))^m(i+hetgen.xcos(j) +  hetgen.zcos(j))))^2))...
%         /((1/4) *dx(i+hetgen.xnos(j))*dz(i+hetgen.znos(j)));
%         end
%     elseif h(i)<0
%      
%             k(i) = sqrt(S(i)) *  ((1-((1- (S(i)^(1/m(i))))^m(i)))^2);
%     else
%         k(i) = 1;
%     end
% 
% end

end