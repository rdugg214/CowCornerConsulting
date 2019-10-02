function [flux_n, flux_w] = CalcFlux(h_p, h_next, k, K, dx)
    flux_w = -k * K_w * ((h_p - h_next)/dx);
    
end