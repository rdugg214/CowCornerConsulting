 function [F] = FVM_index(h,dt, t, t_old,h_old, params, index)
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
t_on_CSG = params{17};
t_on_PUMP = params{18};
S_new = CalcS(h, alpha, n, m);
S_old = CalcS(h_old, alpha, n, m);


k_new = Calck(h, S_new, m);
k_old = Calck(h_old, S_old, m);

psi_new = CalcPsi(h, S_new, psi_res, psi_sat);
psi_old = CalcPsi(h_old, S_old, psi_res, psi_sat);
[Q,Kzz]=  Calc_Q(h,x,z,dt,psi_new,psi_sat,t,t_on_PUMP,Kzz,DELTAX,DELTAZ); %To be updated
Q_old =  Calc_Q(h_old,x,z,dt,psi_old,psi_sat,t,t_on_PUMP,Kzz,DELTAX,DELTAZ);
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
H = h+z;
H_old = h_old+z;
F = zeros(size(h));

%% INside 

i = floor((index-1)/Nx) + 1;
j = mod(index - 1, Nx) + 1;
if (index ~= (i-1)*Nx + j)
    fprintf("input index: %d, calculated index: %d \n", index, (i-1)*Nx + j);
    fprintf("i: %d, j: %d \n", i, j);
    error("This is where is should fucking stop")
end

if i == 1 && j == 1
    flux_n(index) = -((k_new(index+Nx) + k_new(index))/2) * Kzz(index) * ( (H(index+Nx) - H(index))/dz(index));
    flux_e(index) = -((k_new(index+1) + k_new(index))/2) * Kxx(index) * ( (H(index+1) - H(index))/dx(index));
    flux_n_old(index) = -((k_old(index+Nx) + k_old(index))/2) * Kzz(index) * ( (H_old(index+Nx) - H_old(index))/dz(index));
    flux_e_old(index) = -((k_old(index+1)+  k_old(index))/2) * Kxx(index) * ( (H_old(index+1) - H_old(index))/dx(index));
elseif i == 1 && j >=2 && j<= Nx-1
    flux_n(index) = -((k_new(index+Nx) + k_new(index))/2) * Kzz(index) * ( (H(index+Nx) - H(index))/dz(index));
    flux_e(index) = -((k_new(index+1) + k_new(index))/2) * Kxx(index) * ( (H(index+1) - H(index))/dx(index));
    flux_n_old(index) = -((k_old(index+Nx) + k_old(index))/2) * Kzz(index) * ( (H_old(index+Nx) - H_old(index))/dz(index));
    flux_e_old(index) = -((k_old(index+1)+  k_old(index))/2) * Kxx(index) * ( (H_old(index+1) - H_old(index))/dx(index));
elseif  i == 1 && j == Nx
    flux_e(index) = Calc_CSGBound(z(index),h(index),Kxx(index),t,t_on_CSG);
    flux_e_old(index)= Calc_CSGBound(z(index),h_old(index),Kxx(index),t,t_on_CSG);
    flux_n(index) = -((k_new(index+Nx) + k_new(index))/2) * Kzz(index) * ( (H(index+Nx) - H(index))/dz(index));
    flux_w(index) = -((k_new(index)+ k_new(index-1))/2) * Kxx(index) * ( (H(index-1) - H(index))/dx(index-1));
    flux_n_old(index) = -((k_old(index+Nx) + k_old(index))/2) * Kzz(index) * ( (H_old(index+Nx) - H_old(index))/dz(index));
    flux_w_old(index) = -((k_old(index) + k_old(index-1))/2) * Kxx(index) * ( (H_old(index-1) - H_old(index))/dx(index-1));
elseif i >= 2 &&i <= Nz-1 && j ==1
    if 80 <= z(index) && z(index) <= 100
        flux_w(index) = Calc_RiverBound(z(index),h(index));
        flux_w_old(index) = Calc_RiverBound(z(index),h_old(index));
    end
    flux_n(index) = -((k_new(index+Nx) + k_new(index))/2) * Kzz(index) * ( (H(index+Nx) - H(index))/dz(index));
    flux_s(index) = -((k_new(index) + k_new(index-Nx))/2) * Kzz(index) * ( (H(index-Nx) - H(index))/dz(index-Nx));
    flux_e(index) = -((k_new(index+1) + k_new(index))/2) * Kxx(index) * ( (H(index+1) - H(index))/dx(index));
    flux_n_old(index) = -((k_old(index+Nx) + k_old(index))/2) * Kzz(index) * ( (H_old(index+Nx) - H_old(index))/dz(index));
    flux_s_old(index) = -((k_old(index) + k_old(index-Nx))/2) * Kzz(index) * ( (H_old(index-Nx) - H_old(index))/dz(index-Nx));
    flux_e_old(index) = -((k_old(index+1)+  k_old(index))/2) * Kxx(index) * ( (H_old(index+1) - H_old(index))/dx(index));
elseif i >=2 && i<=Nz-1 && j==Nx
    if 0 < z(index) && z(index) <= 40
        flux_e(index) = Calc_CSGBound(z(index),h(index),Kxx(index),t,t_on_CSG);
        flux_e_old(index) = Calc_CSGBound(z(index),h_old(index),Kxx(index),t,t_on_CSG);
    end
    flux_n(index) = -((k_new(index+Nx) + k_new(index))/2) * Kzz(index) * ( (H(index+Nx) - H(index))/dz(index));
    flux_w(index) = -((k_new(index)+ k_new(index-1))/2) * Kxx(index) * ( (H(index-1) - H(index))/dx(index-1));
    flux_s(index) = -((k_new(index) + k_new(index-Nx))/2) * Kzz(index) * ( (H(index-Nx) - H(index))/dz(index-Nx));
    flux_n_old(index) = -((k_old(index+Nx) + k_old(index))/2) * Kzz(index) * ( (H_old(index+Nx) - H_old(index))/dz(index));
    flux_w_old(index) = -((k_old(index) + k_old(index-1))/2) * Kxx(index) * ( (H_old(index-1) - H_old(index))/dx(index-1));
    flux_s_old(index) = -((k_old(index) + k_old(index-Nx))/2) * Kzz(index) * ( (H_old(index-Nx) - H_old(index))/dz(index-Nx));
elseif i == Nz && j == 1
    flux_w(index) = Calc_RiverBound(z(index),h(index));
    flux_w_old(index) = Calc_RiverBound(z(index),h_old(index));
    flux_n(index) = Calc_RainBound(t,h(index),DELTAX(index));
    flux_n_old(index) = Calc_RainBound(t_old,h_old(index),DELTAX(index));
    flux_s(index) = -((k_new(index) + k_new(index-Nx))/2) * Kzz(index) * ( (H(index-Nx) - H(index))/dz(index-Nx));
    flux_e(index) = -((k_new(index+1) + k_new(index))/2) * Kxx(index) * ( (H(index+1) - H(index))/dx(index));
    flux_s_old(index) = -((k_old(index) + k_old(index-Nx))/2) * Kzz(index) * ( (H_old(index-Nx) - H_old(index))/dz(index-Nx));
    flux_e_old(index) = -((k_old(index+1)+  k_old(index))/2) * Kxx(index) * ( (H_old(index+1) - H_old(index))/dx(index));
elseif i == Nz && j >=2 && j<= Nx-1
    flux_n(index) = Calc_RainBound(t,h(index),DELTAX(index));
    flux_n_old(index) = Calc_RainBound(t_old,h_old(index),DELTAX(index));
    flux_s(index) = -((k_new(index) + k_new(index-Nx))/2) * Kzz(index) * ( (H(index-Nx) - H(index))/dz(index-Nx));
    flux_s_old(index) = -((k_old(index) + k_old(index-Nx))/2) * Kzz(index) * ( (H_old(index-Nx) - H_old(index))/dz(index-Nx));
    flux_w(index) = -((k_new(index)+ k_new(index-1))/2) * Kxx(index) * ( (H(index-1) - H(index))/dx(index-1));
    flux_e(index) = -((k_new(index+1) + k_new(index))/2) * Kxx(index) * ( (H(index+1) - H(index))/dx(index));
    flux_w_old(index) = -((k_old(index) + k_old(index-1))/2) * Kxx(index) * ( (H_old(index-1) - H_old(index))/dx(index-1));
    flux_e_old(index) = -((k_old(index+1)+  k_old(index))/2) * Kxx(index) * ( (H_old(index+1) - H_old(index))/dx(index));
elseif i == Nz && j == Nx
    flux_n(index) = Calc_RainBound(t,h(index),DELTAX(index));
    flux_n_old(index) = Calc_RainBound(t_old,h_old(index),DELTAX(index));
    flux_w(index) = -((k_new(index)+ k_new(index-1))/2) * Kxx(index) * ( (H(index-1) - H(index))/dx(index-1));
    flux_s(index) = -((k_new(index) + k_new(index-Nx))/2) * Kzz(index) * ( (H(index-Nx) - H(index))/dz(index-Nx));
    flux_w_old(index) = -((k_old(index) + k_old(index-1))/2) * Kxx(index) * ( (H_old(index-1) - H_old(index))/dx(index-1));
    flux_s_old(index) = -((k_old(index) + k_old(index-Nx))/2) * Kzz(index) * ( (H_old(index-Nx) - H_old(index))/dz(index-Nx));
    
else
    flux_n(index) = -((k_new(index+Nx) + k_new(index))/2) * Kzz(index) * ( (H(index+Nx) - H(index))/dz(index));
    flux_w(index) = -((k_new(index)+ k_new(index-1))/2) * Kxx(index) * ( (H(index-1) - H(index))/dx(index-1));
    flux_s(index) = -((k_new(index) + k_new(index-Nx))/2) * Kzz(index) * ( (H(index-Nx) - H(index))/dz(index-Nx));
    flux_e(index) = -((k_new(index+1) + k_new(index))/2) * Kxx(index) * ( (H(index+1) - H(index))/dx(index));
    
    flux_n_old(index) = -((k_old(index+Nx) + k_old(index))/2) * Kzz(index) * ( (H_old(index+Nx) - H_old(index))/dz(index));
    flux_w_old(index) = -((k_old(index) + k_old(index-1))/2) * Kxx(index) * ( (H_old(index-1) - H_old(index))/dx(index-1));
    flux_s_old(index) = -((k_old(index) + k_old(index-Nx))/2) * Kzz(index) * ( (H_old(index-Nx) - H_old(index))/dz(index-Nx));
    flux_e_old(index) = -((k_old(index+1)+  k_old(index))/2) * Kxx(index) * ( (H_old(index+1) - H_old(index))/dx(index));
end
F = psi_new(index) - psi_old(index) + ...
    theta * dt * (((1/DELTAX(index)) * (flux_e(index) + ...
    flux_w(index))) + ((1/DELTAZ(index)) * (flux_n(index) + ...
    flux_s(index))) - Q(index) )  + ...
    (1 - theta) * dt * ( (1/DELTAX(index) * (flux_e_old(index) + ...
    flux_w_old(index))) + ((1/DELTAZ(index)) * (flux_n_old(index) + ...
    flux_s_old(index))) - Q_old(index));

% F = F(index);
end