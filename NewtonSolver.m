function h = NewtonSolver(h, h_old, t, t_old, dt, fparms)
    MaxIters = 10;
    func = @(h, h_old) FVM(h, dt, t, t_old, h_old, fparms);
    F = func(h, h_old);
    k = 1;
    J = JacobianFD(func, h, h_old);
    while k < MaxIters
        if (CalculateJacobian(k))
            J = JacobianFD(func, h, h_old);
        end
        h_old = h;
        F = F';
        spy(J);
%         disp(J(:, 100))
        disp(det(J));
        dh = J\(-F);
%         dh = gmres(J,-F);
        disp(size(dh));
        h = h + dh';
        F = func(h, h_old);
        k = k + 1;
    end
end