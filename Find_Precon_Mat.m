function M = Find_Precon_Mat(A, type, omega)
    if type == "Jacobi"
        M = diag(diag(A));
    elseif type == "Gauss-Seidel"
        M = tril(A);
    elseif type == "SOR"
        M = tril(A) - (omega - 1)/omega * diag(diag(A));
    else 
        setup.type = 'nofill';
        setup.milu = 'off';
        setup.droptol = 1e-6;
        M = ilu(A, setup);
    end
end