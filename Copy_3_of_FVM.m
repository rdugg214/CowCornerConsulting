 function [F,h_gain] = FVM(h,dt, t, t_old,h_old, params)
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
simple = params{19};
Pr = params{20};
hetgen = params{21};
prediction_data = params{22};
S_new = CalcS(h, alpha, n, m);
S_old = CalcS(h_old, alpha, n, m);


k_new = Calck(h, S_new, m,x,z,dx,dz,hetgen);
k_old = Calck(h_old, S_old, m,x,z,dx,dz,hetgen);

psi_new = CalcPsi(h, S_new, psi_res, psi_sat,x,z,dx,dz,hetgen);
psi_old = CalcPsi(h_old, S_old, psi_res, psi_sat,x,z,dx,dz,hetgen);

[Q,Kzz]=  Calc_Q(h,x,z,dt,psi_new,psi_sat,t,t_on_PUMP,Kzz,DELTAX,DELTAZ,simple,Pr); %To be updated
Q_old =  Calc_Q(h_old,x,z,dt,psi_old,psi_sat,t,t_on_PUMP,Kzz,DELTAX,DELTAZ,simple,Pr);
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
for i = 1:Nz
    for j = 1:Nx
        index = (i-1)*Nx + j;
        if i == 1 && j == 1
        flux_n(index) = -k_face(k_new,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ( (H(index+Nx) - H(index))/dz(index));
        flux_e(index) = -k_face(k_new,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ( (H(index+1) - H(index))/dx(index));
        flux_n_old(index) = -k_face(k_old,dz,index,Nx)*Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ( (H_old(index+Nx) - H_old(index))/dz(index));
        flux_e_old(index) = -k_face(k_old,dx,index,1)*Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ( (H_old(index+1) - H_old(index))/dx(index));
        elseif i == 1 && j >=2 && j<= Nx-1
            flux_n(index) = -k_face(k_new,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ( (H(index+Nx) - H(index))/dz(index));
        flux_e(index) = -k_face(k_new,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ( (H(index+1) - H(index))/dx(index));
        flux_n_old(index) = -k_face(k_old,dz,index,Nx)*Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ( (H_old(index+Nx) - H_old(index))/dz(index));
        flux_e_old(index) = -k_face(k_old,dx,index,1)*Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ( (H_old(index+1) - H_old(index))/dx(index));
        elseif  i == 1 && j == Nx
        flux_e(index) = Calc_CSGBound(z(index),h(index),Kxx(index),t,t_on_CSG);
        flux_e_old(index)= Calc_CSGBound(z(index),h_old(index),Kxx(index),t,t_on_CSG);
        flux_n(index) = -k_face(k_new,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ( (H(index+Nx) - H(index))/dz(index));
        flux_w(index) = -k_face(k_new,dx,index,-1) *Kxz_face(Kxx,-1,index,hetgen.boundary(index)) * ( (H(index-1) - H(index))/dx(index-1));
        flux_n_old(index) = -k_face(k_old,dz,index,Nx)*Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ( (H_old(index+Nx) - H_old(index))/dz(index));
        flux_w_old(index) = -k_face(k_old,dx,index,-1)*Kxz_face(Kxx,-1,index,hetgen.boundary(index)) * ( (H_old(index-1) - H_old(index))/dx(index-1));
        elseif i >= 2 &&i <= Nz-1 && j ==1
         if 80 <= z(index) && z(index) <= 100
            flux_w(index) = Calc_RiverBound(z(index),h(index),simple);
            flux_w_old(index) = Calc_RiverBound(z(index),h_old(index),simple);
        end
        flux_n(index) = -k_face(k_new,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ( (H(index+Nx) - H(index))/dz(index));
        flux_s(index) = - k_face(k_new,dz,index,-Nx)*Kxz_face(Kzz,-Nx,index,hetgen.boundary(index)) * ( (H(index-Nx) - H(index))/dz(index-Nx));
        flux_e(index) = -k_face(k_new,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ( (H(index+1) - H(index))/dx(index));
        flux_n_old(index) = -k_face(k_old,dz,index,Nx)*Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ( (H_old(index+Nx) - H_old(index))/dz(index));
        flux_s_old(index) = -k_face(k_old,dz,index,-Nx)*Kxz_face(Kzz,-Nx,index,hetgen.boundary(index)) * ( (H_old(index-Nx) - H_old(index))/dz(index-Nx));
        flux_e_old(index) = -k_face(k_old,dx,index,1)*Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ( (H_old(index+1) - H_old(index))/dx(index));
        elseif i >=2 && i<=Nz-1 && j==Nx
             if 0 < z(index) && z(index) <= 40
            flux_e(index) = Calc_CSGBound(z(index),h(index),Kxx(index),t,t_on_CSG);
            flux_e_old(index) = Calc_CSGBound(z(index),h_old(index),Kxx(index),t,t_on_CSG);
            end
        flux_n(index) = -k_face(k_new,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ( (H(index+Nx) - H(index))/dz(index));
        flux_w(index) = -k_face(k_new,dx,index,-1) *Kxz_face(Kxx,-1,index,hetgen.boundary(index)) * ( (H(index-1) - H(index))/dx(index-1));
        flux_s(index) = - k_face(k_new,dz,index,-Nx)*Kxz_face(Kzz,-Nx,index,hetgen.boundary(index)) * ( (H(index-Nx) - H(index))/dz(index-Nx));
        flux_n_old(index) = -k_face(k_old,dz,index,Nx)*Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ( (H_old(index+Nx) - H_old(index))/dz(index));
        flux_w_old(index) = -k_face(k_old,dx,index,-1)*Kxz_face(Kxx,-1,index,hetgen.boundary(index)) * ( (H_old(index-1) - H_old(index))/dx(index-1));
        flux_s_old(index) = -k_face(k_old,dz,index,-Nx)*Kxz_face(Kzz,-Nx,index,hetgen.boundary(index)) * ( (H_old(index-Nx) - H_old(index))/dz(index-Nx));
        elseif i == Nz && j == 1
             flux_w(index) = Calc_RiverBound(z(index),h(index),simple);
        flux_w_old(index) = Calc_RiverBound(z(index),h_old(index),simple);
        flux_n(index) = Calc_RainBound(t,h(index),DELTAX(index),simple,prediction_data);
        flux_n_old(index) = Calc_RainBound(t_old,h_old(index),DELTAX(index),simple,prediction_data);
        flux_s(index) = - k_face(k_new,dz,index,-Nx)*Kxz_face(Kzz,-Nx,index,hetgen.boundary(index)) * ( (H(index-Nx) - H(index))/dz(index-Nx));
        flux_e(index) = -k_face(k_new,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ( (H(index+1) - H(index))/dx(index));
        flux_s_old(index) = -k_face(k_old,dz,index,-Nx)*Kxz_face(Kzz,-Nx,index,hetgen.boundary(index)) * ( (H_old(index-Nx) - H_old(index))/dz(index-Nx));
        flux_e_old(index) = -k_face(k_old,dx,index,1)*Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ( (H_old(index+1) - H_old(index))/dx(index));
        elseif i == Nz && j >=2 && j<= Nx-1
            flux_n(index) = Calc_RainBound(t,h(index),DELTAX(index),simple,prediction_data);
        flux_n_old(index) = Calc_RainBound(t_old,h_old(index),DELTAX(index),simple,prediction_data);
        flux_s(index) = - k_face(k_new,dz,index,-Nx)*Kxz_face(Kzz,-Nx,index,hetgen.boundary(index)) * ( (H(index-Nx) - H(index))/dz(index-Nx));
        flux_s_old(index) = -k_face(k_old,dz,index,-Nx)*Kxz_face(Kzz,-Nx,index,hetgen.boundary(index)) * ( (H_old(index-Nx) - H_old(index))/dz(index-Nx));
        flux_w(index) = -k_face(k_new,dx,index,-1) *Kxz_face(Kxx,-1,index,hetgen.boundary(index)) * ( (H(index-1) - H(index))/dx(index-1));
        flux_e(index) = -k_face(k_new,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ( (H(index+1) - H(index))/dx(index));
        flux_w_old(index) = -k_face(k_old,dx,index,-1)*Kxz_face(Kxx,-1,index,hetgen.boundary(index)) * ( (H_old(index-1) - H_old(index))/dx(index-1));
        flux_e_old(index) = -k_face(k_old,dx,index,1)*Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ( (H_old(index+1) - H_old(index))/dx(index));
        elseif i == Nz && j == Nx
         flux_n(index) = Calc_RainBound(t,h(index),DELTAX(index),simple,prediction_data);
        flux_n_old(index) = Calc_RainBound(t_old,h_old(index),DELTAX(index),simple,prediction_data);
        flux_w(index) = -k_face(k_new,dx,index,-1) *Kxz_face(Kxx,-1,index,hetgen.boundary(index)) * ( (H(index-1) - H(index))/dx(index-1));
        flux_s(index) = - k_face(k_new,dz,index,-Nx)*Kxz_face(Kzz,-Nx,index,hetgen.boundary(index)) * ( (H(index-Nx) - H(index))/dz(index-Nx));
        flux_w_old(index) = -k_face(k_old,dx,index,-1)*Kxz_face(Kxx,-1,index,hetgen.boundary(index)) * ( (H_old(index-1) - H_old(index))/dx(index-1));
        flux_s_old(index) = -k_face(k_old,dz,index,-Nx)*Kxz_face(Kzz,-Nx,index,hetgen.boundary(index)) * ( (H_old(index-Nx) - H_old(index))/dz(index-Nx));
   
        else
        flux_n(index) = -k_face(k_new,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ( (H(index+Nx) - H(index))/dz(index));
        flux_w(index) = -k_face(k_new,dx,index,-1) *Kxz_face(Kxx,-1,index,hetgen.boundary(index)) * ( (H(index-1) - H(index))/dx(index-1));
        flux_s(index) = - k_face(k_new,dz,index,-Nx)*Kxz_face(Kzz,-Nx,index,hetgen.boundary(index)) * ( (H(index-Nx) - H(index))/dz(index-Nx));
        flux_e(index) = -k_face(k_new,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ( (H(index+1) - H(index))/dx(index));
       
        flux_n_old(index) = -k_face(k_old,dz,index,Nx)*Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ( (H_old(index+Nx) - H_old(index))/dz(index));
        flux_w_old(index) = -k_face(k_old,dx,index,-1)*Kxz_face(Kxx,-1,index,hetgen.boundary(index)) * ( (H_old(index-1) - H_old(index))/dx(index-1));
        flux_s_old(index) = -k_face(k_old,dz,index,-Nx)*Kxz_face(Kzz,-Nx,index,hetgen.boundary(index)) * ( (H_old(index-Nx) - H_old(index))/dz(index-Nx));
        flux_e_old(index) = -k_face(k_old,dx,index,1)*Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ( (H_old(index+1) - H_old(index))/dx(index));
        end
        F(index) = psi_new(index) - psi_old(index) + ...
            theta * dt * (((1/DELTAX(index)) * (flux_e(index) + ...
            flux_w(index))) + ((1/DELTAZ(index)) * (flux_n(index) + ...
            flux_s(index))) - Q(index) )  + ...
            (1 - theta) * dt * ( (1/DELTAX(index) * (flux_e_old(index) + ...
            flux_w_old(index))) + ((1/DELTAZ(index)) * (flux_n_old(index) + ...
            flux_s_old(index))) - Q_old(index));
    end
end
h_gain.n = flux_n ;
h_gain.e = flux_e ;
h_gain.s = flux_s ;
h_gain.w = flux_w ;
h_gain.Q = Q-Q_old;
% figure(2)

% subplot(1,2,1)
% plot3(x,z,flux_n,'.')
% subplot(1,2,2)
% plot3(x,z,flux_s,'.')
end