function J = JacobianFD(y, h, t, t_old)
%JACOBIAN Returns the Finite Difference Jacobian of F(x) evaluated at x

N = size(y,1);
% J = sparse(N,N);
J = zeros(N, N);
Fx = NewstonSolver(y, x, h, t, t_old);
h_shift = sqrt(eps);
if ~isempty(y)
    h_shift = h_shift*norm(y, 2);
end

for i = 1:N
    xi = y;
    xi(i) = xi(i) + h_shift;
    J(:, i) = (NewstonSolver(xi, x, h, t, t_old) - Fx)/h_shift;
end

end