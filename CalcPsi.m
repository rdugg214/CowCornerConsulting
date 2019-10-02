function psi = CalcPsi(h, S, psi_res, psi_sat)
    psi = ones(size(h));
    for i = 1:length(h)
        if h < 0
            psi(i) = psi_res(i) + S(i)*(psi_sat(i) - psi_res(i));
            
        else
            psi(i) = psi_sat(i);
        end
        
    end
end
