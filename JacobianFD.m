function J = JacobianFD(Ffunc, h, h_old)
%JACOBIAN Returns the Finite Difference Jacobian of F(x) evaluated at x

% N = size(h, 2);
% J = zeros(N, N);
% Fx = Ffunc(h, h_old);
% h_shift = sqrt(eps);
% if ~isempty(h)
%     h_shift = h_shift*norm(h, 2);
% end
% 
% for i = 1:N
%     xi = h;
%     xi(i) = xi(i) + h_shift;
% %     disp(i)
%     J(:, i) = (Ffunc(xi, h_old) - Fx)/h_shift;
% end

x =h;
N = size(h,2);
J = zeros(N,N);
Fx = Ffunc(h, h_old);
if norm(x,2) == 0
    h_shift = sqrt(eps);
else
    h_shift = sqrt(eps)*norm(x,2);
end

e = zeros(N,1);
for j = 1:N
    e(j) = 1.0;
    J(:,j) = (Ffunc(x'+h_shift*e,h_old) - Fx')/h_shift;
    e(j) = 0.0;
end
end