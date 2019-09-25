% Simiple version of the problem

% Define static variables
dt = 1; endtime = 100;

iter_num = 0;
J = JacobianFD(psi, h, t, t_old);
for i = dt:dt:endtime
    iter_num = iter_num + 1;
    if (CalculateJacobian(iter_num))
        J = JacobianFD(psi, h, i, t_old);
    end
    dpsi = J\(-F);
    psi = psi + dpsi;
    NewtonSolver(psi, h, t, t_old)
    t_old = i;
end