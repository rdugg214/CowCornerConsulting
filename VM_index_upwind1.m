 function [F] = FVM_index(h,dt, t, t_old,h_old, params, index, k_old, psi_old, Q_old, k_new, psi_new, Q, Kzz)
% Parameters IN  , k_new, k_old, psi_new, psi_old, Q, Q_old, Kzz, H, H_old
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
% Kzz = params{12};       %vector size h
dx = params{13};        %vector size h
dz = params{14};        %vector size h
DELTAX =  params{15};   %vector size h
DELTAZ = params{16};    %vector size h
t_on_CSG = params{17};
t_on_PUMP = params{18};
simple = params{19};
Pr = params{20};
hetgen = params{21};
sigmax = 1;
sigmaz = 1;
theta = 0.5;

%% Flux Terms INSIDE
flux_n = 0;
flux_w = 0;
flux_n_old = 0;
flux_w_old = 0;
flux_e = 0;
flux_s = 0;
flux_e_old = 0;
flux_s_old = 0;
H = h+z;
H_old = h_old+z;
% F = zeros(size(h));

%% INside 

i = floor((index-1)/Nx) + 1;
j = mod(index - 1, Nx) + 1;
if (index ~= (i-1)*Nx + j)
    fprintf("input index: %d, calculated index: %d \n", index, (i-1)*Nx + j);
    fprintf("i: %d, j: %d \n", i, j);
    error("This is where is should fucking stop")
end
if i>=2 && i<= Nz-1 && j >=2 && j<= Nx-1
    flux_n = - k_face(k_new,dz,index,Nx) * Kzz(index) * ( (H(index) + (sigmaz/2)*(H(index+Nx) - H(index)))/dz(index));
    flux_w = - k_face(k_new,dx,index,-1) * Kxx(index) * ( (H(index-1) - (sigmax/2)*(H(index-1) - H(index)))/dx(index-1));
    flux_s = - k_face(k_new,dz,index,-Nx)* Kzz(index) * ( (H(index-Nx) - (sigmaz/2)*(H(index-Nx) - H(index)))/dz(index-Nx));
    flux_e = - k_face(k_new,dx,index,1)* Kxx(index) * ( (H(index) + (sigmax/2)*(H(index+1) - H(index)))/dx(index));
    
    flux_n_old = - k_face(k_old,dz,index,Nx) * Kzz(index) * ( (H_old(index+Nx) - H_old(index))/dz(index));
    flux_w_old = - k_face(k_old,dx,index,-1) * Kxx(index) * ( (H_old(index-1) - H_old(index))/dx(index-1));
    flux_s_old = - k_face(k_old,dz,index,-Nx) * Kzz(index) * ( (H_old(index-Nx) - H_old(index))/dz(index-Nx));
    flux_e_old = - k_face(k_old,dx,index,1) * Kxx(index) * ( (H_old(index+1) - H_old(index))/dx(index));
elseif i == 1 && j == 1
    flux_n = - k_face(k_new,dz,index,Nx) * Kzz(index) * ( (H(index) + (sigmaz/2)*(H(index+Nx) - H(index)))/dz(index));
    flux_e = - k_face(k_new,dx,index,1)* Kxx(index) * ( (H(index) + (sigmax/2)*(H(index+1) - H(index)))/dx(index));
    flux_n_old = - k_face(k_old,dz,index,Nx) * Kzz(index) * ( (H_old(index+Nx) - H_old(index))/dz(index));
    flux_e_old = - k_face(k_old,dx,index,1) * Kxx(index) * ( (H_old(index+1) - H_old(index))/dx(index));
elseif i == 1 && j >=2 && j<= Nx-1
    flux_n = - k_face(k_new,dz,index,Nx) * Kzz(index) * ( (H(index) + (sigmaz/2)*(H(index+Nx) - H(index)))/dz(index));
    flux_e = - k_face(k_new,dx,index,1)* Kxx(index) * ( (H(index) + (sigmax/2)*(H(index+1) - H(index)))/dx(index));
    flux_n_old = - k_face(k_old,dz,index,Nx) * Kzz(index) * ( (H_old(index+Nx) - H_old(index))/dz(index));
    flux_e_old = - k_face(k_old,dx,index,1) * Kxx(index) * ( (H_old(index+1) - H_old(index))/dx(index));
elseif  i == 1 && j == Nx
    flux_e = Calc_CSGBound(z(index),h(index),Kxx(index),t,t_on_CSG);
    flux_e_old= Calc_CSGBound(z(index),h_old(index),Kxx(index),t,t_on_CSG);
    flux_n = - k_face(k_new,dz,index,Nx) * Kzz(index) * ( (H(index) + (sigmaz/2)*(H(index+Nx) - H(index)))/dz(index));
    flux_w = - k_face(k_new,dx,index,-1) * Kxx(index) * ( (H(index-1) - (sigmax/2)*(H(index-1) - H(index)))/dx(index-1));
    flux_n_old = - k_face(k_old,dz,index,Nx) * Kzz(index) * ( (H_old(index+Nx) - H_old(index))/dz(index));
    flux_w_old = - k_face(k_old,dx,index,-1) * Kxx(index) * ( (H_old(index-1) - H_old(index))/dx(index-1));
elseif i >= 2 &&i <= Nz-1 && j ==1
    if 80 <= z(index) && z(index) <= 100
        flux_w = Calc_RiverBound(z(index),h(index),simple);
        flux_w_old = Calc_RiverBound(z(index),h_old(index),simple);
    end
    flux_n = - k_face(k_new,dz,index,Nx) * Kzz(index) * ( (H(index) + (sigmaz/2)*(H(index+Nx) - H(index)))/dz(index));
    flux_s = - k_face(k_new,dz,index,-Nx)* Kzz(index) * ( (H(index-Nx) - (sigmaz/2)*(H(index-Nx) - H(index)))/dz(index-Nx));
    flux_e = - k_face(k_new,dx,index,1)* Kxx(index) * ( (H(index) + (sigmax/2)*(H(index+1) - H(index)))/dx(index));
    flux_n_old = - k_face(k_old,dz,index,Nx) * Kzz(index) * ( (H_old(index+Nx) - H_old(index))/dz(index));
    flux_s_old = - k_face(k_old,dz,index,-Nx) * Kzz(index) * ( (H_old(index-Nx) - H_old(index))/dz(index-Nx));
    flux_e_old = - k_face(k_old,dx,index,1) * Kxx(index) * ( (H_old(index+1) - H_old(index))/dx(index));
elseif i >=2 && i<=Nz-1 && j==Nx
    if 0 < z(index) && z(index) <= 40
        flux_e = Calc_CSGBound(z(index),h(index),Kxx(index),t,t_on_CSG);
        flux_e_old = Calc_CSGBound(z(index),h_old(index),Kxx(index),t,t_on_CSG);
    end
    flux_n = - k_face(k_new,dz,index,Nx) * Kzz(index) * ( (H(index) + (sigmaz/2)*(H(index+Nx) - H(index)))/dz(index));
    flux_w = - k_face(k_new,dx,index,-1) * Kxx(index) * ( (H(index-1) - (sigmax/2)*(H(index-1) - H(index)))/dx(index-1));
    flux_s = - k_face(k_new,dz,index,-Nx)* Kzz(index) * ( (H(index-Nx) - (sigmaz/2)*(H(index-Nx) - H(index)))/dz(index-Nx));
    flux_n_old = - k_face(k_old,dz,index,Nx) * Kzz(index) * ( (H_old(index+Nx) - H_old(index))/dz(index));
    flux_w_old = - k_face(k_old,dx,index,-1) * Kxx(index) * ( (H_old(index-1) - H_old(index))/dx(index-1));
    flux_s_old = - k_face(k_old,dz,index,-Nx) * Kzz(index) * ( (H_old(index-Nx) - H_old(index))/dz(index-Nx));
elseif i == Nz && j == 1
    flux_w = Calc_RiverBound(z(index),h(index),simple);
    flux_w_old = Calc_RiverBound(z(index),h_old(index),simple);
    flux_n = Calc_RainBound(t,h(index),DELTAX(index));
    flux_n_old = Calc_RainBound(t_old,h_old(index),DELTAX(index));
    flux_s = - k_face(k_new,dz,index,-Nx)* Kzz(index) * ( (H(index-Nx) - (sigmaz/2)*(H(index-Nx) - H(index)))/dz(index-Nx));
    flux_e = - k_face(k_new,dx,index,1)* Kxx(index) * ( (H(index) + (sigmax/2)*(H(index+1) - H(index)))/dx(index));
    flux_s_old = - k_face(k_old,dz,index,-Nx) * Kzz(index) * ( (H_old(index-Nx) - H_old(index))/dz(index-Nx));
    flux_e_old = - k_face(k_old,dx,index,1) * Kxx(index) * ( (H_old(index+1) - H_old(index))/dx(index));
elseif i == Nz && j >=2 && j<= Nx-1
    flux_n = Calc_RainBound(t,h(index),DELTAX(index));
    flux_n_old = Calc_RainBound(t_old,h_old(index),DELTAX(index));
    flux_s = - k_face(k_new,dz,index,-Nx)* Kzz(index) * ( (H(index-Nx) - (sigmaz/2)*(H(index-Nx) - H(index)))/dz(index-Nx));
    flux_s_old = - k_face(k_old,dz,index,-Nx) * Kzz(index) * ( (H_old(index-Nx) - H_old(index))/dz(index-Nx));
    flux_w = - k_face(k_new,dx,index,-1) * Kxx(index) * ( (H(index-1) - (sigmax/2)*(H(index-1) - H(index)))/dx(index-1));
    flux_e = - k_face(k_new,dx,index,1)* Kxx(index) * ( (H(index) + (sigmax/2)*(H(index+1) - H(index)))/dx(index));
    flux_w_old = - k_face(k_old,dx,index,-1) * Kxx(index) * ( (H_old(index-1) - H_old(index))/dx(index-1));
    flux_e_old = - k_face(k_old,dx,index,1) * Kxx(index) * ( (H_old(index+1) - H_old(index))/dx(index));
elseif i == Nz && j == Nx
    flux_n = Calc_RainBound(t,h(index),DELTAX(index));
    flux_n_old = Calc_RainBound(t_old,h_old(index),DELTAX(index));
    flux_w = - k_face(k_new,dx,index,-1) * Kxx(index) * ( (H(index-1) - (sigmax/2)*(H(index-1) - H(index)))/dx(index-1));
    flux_s = - k_face(k_new,dz,index,-Nx)* Kzz(index) * ( (H(index-Nx) - (sigmaz/2)*(H(index-Nx) - H(index)))/dz(index-Nx));
    flux_w_old = - k_face(k_old,dx,index,-1) * Kxx(index) * ( (H_old(index-1) - H_old(index))/dx(index-1));
    flux_s_old = - k_face(k_old,dz,index,-Nx) * Kzz(index) * ( (H_old(index-Nx) - H_old(index))/dz(index-Nx));
end
F = psi_new(index) - psi_old(index) + ...
    theta * dt * (((1/DELTAX(index)) * (flux_e + ...
    flux_w)) + ((1/DELTAZ(index)) * (flux_n + ...
    flux_s)) - Q(index) )  + ...
    (1 - theta) * dt * ( (1/DELTAX(index) * (flux_e_old + ...
    flux_w_old)) + ((1/DELTAZ(index)) * (flux_n_old + ...
    flux_s_old)) - Q_old(index));
% F = F(index);
end