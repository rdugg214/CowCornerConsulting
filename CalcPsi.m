function psi = CalcPsi(h, S, psi_res, psi_sat)
psi = ones(size(h)).*psi_sat;
locs = h<0;
% psi(locs) = psi_res(locs) + S(locs).*(psi_sat(locs) - psi_res(locs));
% psi(~locs) = psi_sat(~locs);
% for i = 1:length(h)
%     if h(i)<0 
%     psi(i) = psi_res(i)+S(i) *(psi_sat(i) - psi_res(i));
% end
psi(locs) = psi_res(locs) + S(locs).*(psi_sat(locs)-psi_res(locs));
end
