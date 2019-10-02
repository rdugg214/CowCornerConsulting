function F = FVM(h,dt, t, t_old,h_old, params)

% Parameters IN
N = params{1};          %value
Nx = params{2};         %value
Nz = params{3};         %value
alpha = params{4};      %vector size h
n = params{5};          %vector size h
m = params{6};          %vector size h
psi_res = params{7};    %vector size h
psi_sat = params{8};    %vector size h
x = params{9};          %vector size h
z = params{10};         %vector size h
Kxx = params{11};       %vector size h
Kzz = params{12};       %vector size h
dx = params{13};        %vector size h
dz = params{14};        %vector size h
DELTAX =  params{15};   %vector size h
DELTAZ = params{16};    %vector size h
alpha = params{17};
% Sub functions and calculations IN
S_new = CalcS(h, alpha, n, m);
S_old = CalcS(h_old, alpha, n, m);
k_new = Calck(h, S_new, m);
k_old = Calck(h_old, S_old, m);
psi_new = CalcPsi(h, S_new, psi_res, psi_sat);
psi_old = CalcPsi(h_old, S_old, psi_res, psi_sat);
Q = CalcQ(h,x,z); %To be updated
Q_old = CalcQ(h_old,x,z);
theta = 0.5;

%% Flux Terms INSIDE
flux_n = zeros(size(x));
flux_w = zeros(size(x));
flux_n_old = zeros(size(x));
flux_w_old = zeros(size(x));
flux_e = zeros(size(x));
flux_s = zeros(size(x));
flux_e_old = zeros(size(x));
flux_s_old = zeros(size(x));
for i = 2:Nz-1
    for j = 2:Nx-1
        index = (i-1)*Nx + j;
        flux_n(index) = -k_new(index) * Kzz(index) * ((h(index+Nx) - h(index))/dz(index) + 1);
        flux_w(index) = -k_new(index) * Kxx(index) * ((h(index) - h(index-1))/dx(index));
        flux_s(index) = -k_new(index) * Kzz(index) * ((h(index) - h(index-Nx))/dz(index) + 1);
        flux_e(index) = -k_new(index) * Kxx(index) * ((h(index+1) - h(index))/dx(index));
        flux_n_old(index) = -k_old(index) * Kzz(index) * ((h_old(index+Nx) - h_old(index))/dz(index) + 1);
        flux_w_old(index) = -k_old(index) * Kxx(index) * ((h_old(index) - h_old(index-1))/dx(index));
        flux_s_old(index) = -k_new(index) * Kzz(index) * ((h_old(index) - h_old(index-Nx))/dz(index) + 1);
        flux_e_old(index) = -k_new(index) * Kxx(index) * ((h_old(index+1) - h_old(index))/dx(index));
    end
end


%% Region 1
F = zeros(size(h));
for i = 2:Nz-1 %loop by row (left to right, dowtn to top)
    for j = 2:Nx-1
        index = (i-1)*Nx + j;
        
        F(index) = psi_new(index) - psi_old(index) + ...
            theta * dt * (1/DELTAX(index) * (flux_e(index) - ...
            flux_w(index)) + 1/DELTAZ(index) * (flux_n(index) - ...
            flux_s(index)) - Q(index)) + ...
            (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index) - ...
            flux_w_old(index)) + 1/DELTAZ(index) * (flux_n_old(index) - ...
            flux_s_old(index)) - Q_old(index));
    end
end


%% Region 2
for i = Nz
    for j = 1
        index = (i-1)*Nx + j;
        
        fluxw_river = Calc_RiverBound(z(index),h(index));
        fluxw_river_old = Calc_RiverBound(z(index),h_old(index));
        fluxn_rain = Calc_RainBound(t,h(index));
        fluxn_rain_old = Calc_RainBound(t_old,h_old(index));
        
        F(index) = psi_new(index) - psi_old(index) + ...
            theta * dt * (1/DELTAX(index) * (flux_e(index) - ...
            fluxw_river) + 1/DELTAZ(index) * (fluxn_rain - ...
            flux_s(index)) - Q(index)) + ...
            (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index) - ...
            fluxw_river_old) + 1/DELTAZ(index) * (fluxn_rain_old - ...
            flux_s_old(index)) - Q_old(index));
    end
end

%% Region 7
for i = Nz
    for j = Nx
        index = (i-1)*Nx + j;
        
        
        fluxn_rain = Calc_RainBound(t,h(index));
        fluxn_rain_old = Calc_RainBound(t_old,h_old(index));
        
        F(index) = psi_new(index) - psi_old(index) + ...
            theta * dt * (1/DELTAX(index) * ( - ...
            flux_w(index)) + 1/DELTAZ(index) * (fluxn_rain - ...
            flux_s(index)) - Q(index)) + ...
            (1 - theta) * dt * (1/DELTAX(index) * ( - ...
            flux_w_old(index)) + 1/DELTAZ(index) * (fluxn_rain_old - ...
            flux_s_old(index)) - Q_old(index));
    end
end

%% Region 4
for i = 1
    for j = 1
        index = (i-1)*Nx + j;
       
        
        F(index) = psi_new(index) - psi_old(index) + ...
            theta * dt * (1/DELTAX(index) * ( ...
            flux_e(index)) + 1/DELTAZ(index) * ( flux_n(index) ...
            ) - Q(index)) + ...
            (1 - theta) * dt * (1/DELTAX(index) * (  ...
            flux_e_old(index)) + 1/DELTAZ(index) * ( ...
            flux_n_old(index)) - Q_old(index));
    end
end

%% Region 9
for i = 1
    for j = Nx
        index = (i-1)*Nx + j;
         fluxe_csg = Calc_CSGBound(z(index),h(index),Kxx(index));
        fluxe_csg_old = Calc_RiverBound(z(index),h_old(index));
        
        F(index) = psi_new(index) - psi_old(index) + ...
            theta * dt * (1/DELTAX(index) * (fluxe_csg- ...
            flux_e(index)) + 1/DELTAZ(index) * ( flux_n(index) ...
            ) - Q(index)) + ...
            (1 - theta) * dt * (1/DELTAX(index) * (  ...
            flux_e_old(index)) + 1/DELTAZ(index) * (fluxe_csg_old- ...
            flux_n_old(index)) - Q_old(index));
    end
end

%%  Region 3
for i = 2:Nz-1
    for j = 1
        index = (i-1)*Nx + j;
        
        if 80 <= z(index) & z(index) <= 100
        FLUX_W_Flat = Calc_RiverBound(z(index),h(index));
        FLUX_W_Flat_old = Calc_RiverBound(z(index),h_old(index));
        else
        FLUX_W_Flat = 0;
        FLUX_W_Flat_old = 0;
        end
        
        F(index) = psi_new(index) - psi_old(index) + ...
            theta * dt * (1/DELTAX(index) * (flux_e(index) - ...
            FLUX_W_Flat ) + 1/DELTAZ(index) * ( flux_n(index) - ...
            flux_s(index)) - Q(index)) + ...
            (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index) - ...
            FLUX_W_Flat_old ) + 1/DELTAZ(index) * (flux_n_old(index) - ...
            flux_s_old(index)) - Q_old(index));
    end
end
% 
%%  Region 8
for i = 2:Nz-1
    for j = Nx
        index = (i-1)*Nx + j;
        
        if 0 < z(index) & z(index) <= 5
        FLUX_E_Flat = Calc_RiverBound(z(index),h(index));
        FLUX_E_Flat_old = Calc_RiverBound(z(index),h_old(index));
        else
        FLUX_E_Flat = 0;
        FLUX_E_Flat_old = 0;
        end
        
        F(index) = psi_new(index) - psi_old(index) + ...
            theta * dt * (1/DELTAX(index) * (FLUX_E_Flat - ...
            flux_w(index) ) + 1/DELTAZ(index) * (flux_n(index) - ...
            flux_s(index)) - Q(index)) + ...
            (1 - theta) * dt * (1/DELTAX(index) * (FLUX_E_Flat_old - ...
            flux_w_old(index)) + 1/DELTAZ(index) * (flux_n_old(index) - ...
            flux_s_old(index)) - Q_old(index));
    end
end

%% Region 5
for i = Nz
    for j = 2:Nx-1
        index = (i-1)*Nx + j;
        
        
        fluxn_rain = Calc_RainBound(t,h(index));
        fluxn_rain_old = Calc_RainBound(t_old,h_old(index));
        
        F(index) = psi_new(index) - psi_old(index) + ...
            theta * dt * (1/DELTAX(index) * (flux_e(index) - ...
            flux_w(index)) + 1/DELTAZ(index) * (fluxn_rain - ...
            flux_s(index)) - Q(index)) + ...
            (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index) - ...
            flux_w_old(index)) + 1/DELTAZ(index) * (fluxn_rain_old - ...
            flux_s_old(index)) - Q_old(index));
    end
end

%% Region 6
for i = 1
    for j = 2:Nx-1
        index = (i-1)*Nx + j;
        
        F(index) = psi_new(index) - psi_old(index) + ...
            theta * dt * (1/DELTAX(index) * (flux_e(index) - ...
            flux_w(index)) + 1/DELTAZ(index) * (flux_n(index)) - Q(index)) + ...
            (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index) - ...
            flux_w_old(index)) + 1/DELTAZ(index) * (flux_n_old(index)) - Q_old(index));
    end
end

% disp('SIM DONE')

end