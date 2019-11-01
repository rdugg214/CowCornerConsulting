function psi = CalcPsi(h, S, psi_res, psi_sat,x,z,dx,dz,hetgen)
%% Homogenous Psi Calc
psi = ones(size(h)).*psi_sat;
locs = h<0;
psi(locs) = psi_res(locs) + S(locs).*(psi_sat(locs)-psi_res(locs));

%% Heterogeneous
% for i = 1:length(h)
%     if h(i)<0 && hetgen.boundary(i) && (z(i) == 100 || z(i) == 0 || x(i) == 0 || x(i) == 500)
%      psi(i) = psi_res(i)+S(i) *(psi_sat(i) - psi_res(i));
%     elseif h(i)<0 && hetgen.boundary(i)
%          psi(i) = 0;
%         for j = 1:4
%      psi(i) = psi(i) + (psi_res(i+hetgen.xcos(j) +  hetgen.zcos(j))+S(i) ...
%          *(psi_sat(i+hetgen.xcos(j) +  hetgen.zcos(j)) - psi_res(i+hetgen.xcos(j) +  hetgen.zcos(j))))/...
%          ((1/4) *dx(i+hetgen.xnos(j))*dz(i+hetgen.znos(j)));
%         end
%     elseif h(i)<0 
%     psi(i) = psi_res(i)+S(i) *(psi_sat(i) - psi_res(i));
%     else
%     psi(i) = psi_sat(i);
%     end
% end

end
