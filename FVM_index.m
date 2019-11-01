 function [F] = FVM_index(h,dt, t, t_old,h_old, params, index_list, k_old, psi_old, Q_old, k_new, psi_new, Q, Kzz)
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
prediction_data = params{22};
DELCSG = params{23};
rain = params{24};
theta = 0.5;

%% Flux Terms INSIDE
flux_n = zeros(size(x));
flux_n_old = zeros(size(x));
flux_e = zeros(size(x));
flux_e_old = zeros(size(x));
H = h+z;
H_old = h_old+z;

F = zeros(size(h));


%% INside 

num_indexes = length(index_list);
i_list = zeros(size(index_list));
j_list = zeros(size(index_list));

for k = 1:num_indexes
    i_list(k) = floor((index_list(k)-1)/Nx) + 1;
    j_list(k) = mod(index_list(k) - 1, Nx) + 1;
end

for k = 1:num_indexes
    index = index_list(k);
    i = i_list(k);
    j = j_list(k);
    if i>=2 && i<= Nz-1 && j >=2 && j<= Nx-1 % Internal
        flux_n(index) = -k_face(k_new,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H(index+Nx) - H(index))/dz(index));
        flux_e(index) = -k_face(k_new,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H(index+1) - H(index))/dx(index));
        flux_n_old(index) = -k_face(k_old,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H_old(index+Nx) - H_old(index))/dz(index));
        flux_e_old(index) = -k_face(k_old,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H_old(index+1) - H_old(index))/dx(index));
    elseif i == 1 && j == 1 % Bottom Left
        flux_n(index) = -k_face(k_new,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H(index+Nx) - H(index))/dz(index));
        flux_e(index) = -k_face(k_new,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H(index+1) - H(index))/dx(index));
        flux_n_old(index) = -k_face(k_old,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H_old(index+Nx) - H_old(index))/dz(index));
        flux_e_old(index) = -k_face(k_old,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H_old(index+1) - H_old(index))/dx(index));
    elseif i == 1 && j >=2 && j<= Nx-1 % Bottom
        flux_n(index) = -k_face(k_new,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H(index+Nx) - H(index))/dz(index));
        flux_e(index) = -k_face(k_new,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H(index+1) - H(index))/dx(index));
        flux_n_old(index) = -k_face(k_old,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H_old(index+Nx) - H_old(index))/dz(index));
        flux_e_old(index) = -k_face(k_old,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H_old(index+1) - H_old(index))/dx(index));
    elseif  i == 1 && j == Nx % Bottom Right
        flux_e(index) = Calc_CSGBound(z(index),h(index),Kxx(index),t,t_on_CSG,DELCSG);
        flux_e_old(index)= Calc_CSGBound(z(index),h_old(index),Kxx(index),t_old,t_on_CSG,DELCSG);
        flux_n(index) = -k_face(k_new,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H(index+Nx) - H(index))/dz(index));
        flux_n_old(index) = -k_face(k_old,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H_old(index+Nx) - H_old(index))/dz(index));
    elseif i >= 2 &&i <= Nz-1 && j ==1 % Left
        flux_n(index) = -k_face(k_new,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H(index+Nx) - H(index))/dz(index));
        flux_e(index) = -k_face(k_new,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H(index+1) - H(index))/dx(index));
        flux_n_old(index) = -k_face(k_old,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H_old(index+Nx) - H_old(index))/dz(index));
        flux_e_old(index) = -k_face(k_old,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H_old(index+1) - H_old(index))/dx(index));
    elseif i >=2 && i<=Nz-1 && j==Nx % Right
        if 0 < z(index) && z(index) <= 5
            flux_e(index) = Calc_CSGBound(z(index),h(index),Kxx(index),t,t_on_CSG,DELCSG);
            flux_e_old(index) = Calc_CSGBound(z(index),h_old(index),Kxx(index),t_old,t_on_CSG,DELCSG);
        end
        flux_n(index) = -k_face(k_new,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H(index+Nx) - H(index))/dz(index));
        flux_n_old(index) = -k_face(k_old,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H_old(index+Nx) - H_old(index))/dz(index));
    elseif i == Nz && j == 1 % Top Left
        flux_w =  Calc_RiverBound(z(index),h(index),simple);
        flux_w_old =  Calc_RiverBound(z(index),h_old(index),simple);
        flux_n(index) = Calc_RainBound(t,h(index),DELTAX(index),rain,prediction_data);
        flux_n_old(index) = Calc_RainBound(t_old,h_old(index),DELTAX(index),rain,prediction_data);
        flux_e(index) = -k_face(k_new,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H(index+1) - H(index))/dx(index));
        flux_e_old(index) = -k_face(k_old,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H_old(index+1) - H_old(index))/dx(index));
    elseif i == Nz && j >=2 && j <= Nx-1 % Top
        flux_n(index) = Calc_RainBound(t,h(index),DELTAX(index),rain,prediction_data);
        flux_n_old(index) = Calc_RainBound(t_old,h_old(index),DELTAX(index),rain,prediction_data);
        flux_e(index) = -k_face(k_new,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H(index+1) - H(index))/dx(index));
        flux_e_old(index) = -k_face(k_old,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H_old(index+1) - H_old(index))/dx(index));
    elseif i == Nz && j == Nx % Top Right
        flux_n(index) = Calc_RainBound(t,h(index),DELTAX(index),rain,prediction_data);
        flux_n_old(index) = Calc_RainBound(t_old,h_old(index),DELTAX(index),rain,prediction_data);
    end
end

for k = 1:num_indexes
    index = index_list(k);
    i = i_list(k);
    j = j_list(k);
    if i>=2 && i<= Nz-1 && j >=2 && j<= Nx-1 % Internal
        F(index) = psi_new(index) - psi_old(index) + ...
            theta * dt * (1/DELTAX(index) * (flux_e(index) - ...
            flux_e(index-1)) + 1/DELTAZ(index) * (flux_n(index) - ...
            flux_n(index-Nx)) - Q(index)) + ...
            (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index) - ...
            flux_e_old(index-1)) + 1/DELTAZ(index) * (flux_n_old(index) - ...
            flux_n_old(index-Nx)) - Q_old(index));
    elseif i == 1 && j == 1 % Bottom Left
        F(index) = psi_new(index) - psi_old(index) + ...
            theta * dt * (1/DELTAX(index) * (flux_e(index) ...
            ) + 1/DELTAZ(index) * (flux_n(index)) - Q(index)) + ...
            (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index)  ...
            ) + 1/DELTAZ(index) * (flux_n_old(index)) - Q_old(index));
    elseif i == 1 && j >=2 && j<= Nx-1 % Bottom
        F(index) = psi_new(index) - psi_old(index) + ...
            theta * dt * (1/DELTAX(index) * (flux_e(index) - ...
            flux_e(index-1)) + 1/DELTAZ(index) * (flux_n(index)) - Q(index)) + ...
            (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index) - ...
            flux_e_old(index-1)) + 1/DELTAZ(index) * (flux_n_old(index)) - Q_old(index));
    elseif  i == 1 && j == Nx % Bottom Right
        F(index) = psi_new(index) - psi_old(index) + ...
            theta * dt * (1/DELTAX(index) * (flux_e(index) - ...
            flux_e(index-1)) + 1/DELTAZ(index) * (flux_n(index)) - Q(index)) + ...
            (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index) - ...
            flux_e_old(index-1)) + 1/DELTAZ(index) * (flux_n_old(index)) - Q_old(index));
    elseif i >= 2 &&i <= Nz-1 && j ==1 % Left
        if 80 <= z(index) && z(index) <= 100
            flux_w =  Calc_RiverBound(z(index),h(index),simple);
            flux_w_old =  Calc_RiverBound(z(index),h_old(index),simple);
            F(index) = psi_new(index) - psi_old(index) + ...
                theta * dt * (1/DELTAX(index) * (flux_e(index) - ...
                flux_w) + 1/DELTAZ(index) * (flux_n(index) - ...
                flux_n(index-Nx)) - Q(index)) + ...
                (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index) - ...
                flux_w_old ) + 1/DELTAZ(index) * (flux_n_old(index) - ...
                flux_n_old(index-Nx)) - Q_old(index));
        else
            F(index) = psi_new(index) - psi_old(index) + ...
                theta * dt * (1/DELTAX(index) * (flux_e(index) ...
                ) + 1/DELTAZ(index) * (flux_n(index) - ...
                flux_n(index-Nx)) - Q(index)) + ...
                (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index)  ...
                ) + 1/DELTAZ(index) * (flux_n_old(index) - ...
                flux_n_old(index-Nx)) - Q_old(index));
        end
    elseif i >=2 && i<=Nz-1 && j==Nx % Right
        F(index) = psi_new(index) - psi_old(index) + ...
            theta * dt * (1/DELTAX(index) * (flux_e(index) - ...
            flux_e(index-1)) + 1/DELTAZ(index) * (flux_n(index) - ...
            flux_n(index-Nx)) - Q(index)) + ...
            (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index) - ...
            flux_e_old(index-1)) + 1/DELTAZ(index) * (flux_n_old(index) - ...
            flux_n_old(index-Nx)) - Q_old(index));
    elseif i == Nz && j == 1 % Top Left
        F(index) = psi_new(index) - psi_old(index) + ...
            theta * dt * (1/DELTAX(index) * (flux_e(index) - ...
            flux_w) + 1/DELTAZ(index) * (flux_n(index) - ...
            flux_n(index-Nx)) - Q(index)) + ...
            (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index) - ...
            flux_w_old) + 1/DELTAZ(index) * (flux_n_old(index) - ...
            flux_n_old(index-Nx)) - Q_old(index));
    elseif i == Nz && j >=2 && j <= Nx-1 % Top
        F(index) = psi_new(index) - psi_old(index) + ...
            theta * dt * (1/DELTAX(index) * (flux_e(index) - ...
            flux_e(index-1)) + 1/DELTAZ(index) * (flux_n(index) - ...
            flux_n(index-Nx)) - Q(index)) + ...
            (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index) - ...
            flux_e_old(index-1)) + 1/DELTAZ(index) * (flux_n_old(index) - ...
            flux_n_old(index-Nx)) - Q_old(index));
    elseif i == Nz && j == Nx % Top Right
        F(index) = psi_new(index) - psi_old(index) + ...
            theta * dt * (1/DELTAX(index) * ( - ...
            flux_e(index-1)) + 1/DELTAZ(index) * (flux_n(index) - ...
            flux_n(index-Nx)) - Q(index)) + ...
            (1 - theta) * dt * (1/DELTAX(index) * ( - ...
            flux_e_old(index-1)) + 1/DELTAZ(index) * (flux_n_old(index) - ...
            flux_n_old(index-Nx)) - Q_old(index));
    end
end
F = F(index_list(num_indexes));
% for k = 1:num_indexes
%     index = index_list(k);
%     if i_list(k)>=2 && i_list(k)<= Nz-1 && j_list(k) >=2 && j_list(k)<= Nx-1 % Internal
%         flux_n(index) = -k_face(k_new,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H(index+Nx) - H(index))/dz(index));
%         flux_e(index) = -k_face(k_new,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H(index+1) - H(index))/dx(index));
%         flux_n_old(index) = -k_face(k_old,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H_old(index+Nx) - H_old(index))/dz(index));
%         flux_e_old(index) = -k_face(k_old,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H_old(index+1) - H_old(index))/dx(index));
%     end
% end
% 
% for k = 1:num_indexes
%     index = index_list(k);
%     if i_list(k) == 1 && j_list(k) == 1 % Bottom Left
%         flux_n(index) = -k_face(k_new,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H(index+Nx) - H(index))/dz(index));
%         flux_e(index) = -k_face(k_new,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H(index+1) - H(index))/dx(index));
%         flux_n_old(index) = -k_face(k_old,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H_old(index+Nx) - H_old(index))/dz(index));
%         flux_e_old(index) = -k_face(k_old,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H_old(index+1) - H_old(index))/dx(index));
%         
%         F(index) = psi_new(index) - psi_old(index) + ...
%             theta * dt * (1/DELTAX(index) * (flux_e(index) ...
%             ) + 1/DELTAZ(index) * (flux_n(index)) - Q(index)) + ...
%             (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index)  ...
%             ) + 1/DELTAZ(index) * (flux_n_old(index)) - Q_old(index));
%     end
% end
% 
% for k = 1:num_indexes
%     index = index_list(k);
%     if i_list(k) == 1 && j_list(k) >=2 && j_list(k)<= Nx-1 % Bottom
%         flux_n(index) = -k_face(k_new,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H(index+Nx) - H(index))/dz(index));
%         flux_e(index) = -k_face(k_new,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H(index+1) - H(index))/dx(index));
%         flux_n_old(index) = -k_face(k_old,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H_old(index+Nx) - H_old(index))/dz(index));
%         flux_e_old(index) = -k_face(k_old,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H_old(index+1) - H_old(index))/dx(index));
%         
%         F(index) = psi_new(index) - psi_old(index) + ...
%             theta * dt * (1/DELTAX(index) * (flux_e(index) - ...
%             flux_e(index-1)) + 1/DELTAZ(index) * (flux_n(index)) - Q(index)) + ...
%             (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index) - ...
%             flux_e_old(index-1)) + 1/DELTAZ(index) * (flux_n_old(index)) - Q_old(index));
%     end
% end
% 
% for k = 1:num_indexes
%     index = index_list(k);
%     if  i_list(k) == 1 && j_list(k) == Nx % Bottom Right
%         flux_e(index) = Calc_CSGBound(z(index),h(index),Kxx(index),t,t_on_CSG,DELCSG);
%         flux_e_old(index)= Calc_CSGBound(z(index),h_old(index),Kxx(index),t_old,t_on_CSG,DELCSG);
%         
%         flux_n(index) = -k_face(k_new,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H(index+Nx) - H(index))/dz(index));
%         flux_n_old(index) = -k_face(k_old,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H_old(index+Nx) - H_old(index))/dz(index));
%         
%         F(index) = psi_new(index) - psi_old(index) + ...
%             theta * dt * (1/DELTAX(index) * (flux_e(index) - ...
%             flux_e(index-1)) + 1/DELTAZ(index) * (flux_n(index)) - Q(index)) + ...
%             (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index) - ...
%             flux_e_old(index-1)) + 1/DELTAZ(index) * (flux_n_old(index)) - Q_old(index));
%     end
% end
% 
% for k = 1:num_indexes
%     index = index_list(k);
%     if i_list(k) >= 2 &&i_list(k) <= Nz-1 && j_list(k) ==1 % Left
%         flux_n(index) = -k_face(k_new,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H(index+Nx) - H(index))/dz(index));
%         flux_e(index) = -k_face(k_new,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H(index+1) - H(index))/dx(index));
%         flux_n_old(index) = -k_face(k_old,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H_old(index+Nx) - H_old(index))/dz(index));
%         flux_e_old(index) = -k_face(k_old,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H_old(index+1) - H_old(index))/dx(index));
%         
%         if 80 <= z(index) && z(index) <= 100
%             flux_w =  Calc_RiverBound(z(index),h(index),simple);
%             flux_w_old =  Calc_RiverBound(z(index),h_old(index),simple);
%             F(index) = psi_new(index) - psi_old(index) + ...
%                 theta * dt * (1/DELTAX(index) * (flux_e(index) - ...
%                 flux_w) + 1/DELTAZ(index) * (flux_n(index) - ...
%                 flux_n(index-Nx)) - Q(index)) + ...
%                 (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index) - ...
%                 flux_w_old ) + 1/DELTAZ(index) * (flux_n_old(index) - ...
%                 flux_n_old(index-Nx)) - Q_old(index));
%             
%             
%         else
%             F(index) = psi_new(index) - psi_old(index) + ...
%                 theta * dt * (1/DELTAX(index) * (flux_e(index) ...
%                 ) + 1/DELTAZ(index) * (flux_n(index) - ...
%                 flux_n(index-Nx)) - Q(index)) + ...
%                 (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index)  ...
%                 ) + 1/DELTAZ(index) * (flux_n_old(index) - ...
%                 flux_n_old(index-Nx)) - Q_old(index));
%         end
%     end
% end
% 
% for k = 1:num_indexes
%     index = index_list(k);
%     if i_list(k) >=2 && i_list(k)<=Nz-1 && j_list(k)==Nx % Right
%         if 0 < z(index) && z(index) <= 5
%             flux_e(index) = Calc_CSGBound(z(index),h(index),Kxx(index),t,t_on_CSG,DELCSG);
%             flux_e_old(index) = Calc_CSGBound(z(index),h_old(index),Kxx(index),t_old,t_on_CSG,DELCSG);
%         end
%         flux_n(index) = -k_face(k_new,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H(index+Nx) - H(index))/dz(index));
%         flux_n_old(index) = -k_face(k_old,dz,index,Nx) * Kxz_face(Kzz,Nx,index,hetgen.boundary(index)) * ((H_old(index+Nx) - H_old(index))/dz(index));
%         
%         F(index) = psi_new(index) - psi_old(index) + ...
%             theta * dt * (1/DELTAX(index) * (flux_e(index) - ...
%             flux_e(index-1)) + 1/DELTAZ(index) * (flux_n(index) - ...
%             flux_n(index-Nx)) - Q(index)) + ...
%             (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index) - ...
%             flux_e_old(index-1)) + 1/DELTAZ(index) * (flux_n_old(index) - ...
%             flux_n_old(index-Nx)) - Q_old(index));
%     end
% end
% 
% for k = 1:num_indexes
%     index = index_list(k);
%     if i_list(k) == Nz && j_list(k) == 1 % Top Left
%         flux_w =  Calc_RiverBound(z(index),h(index),simple);
%         flux_w_old =  Calc_RiverBound(z(index),h_old(index),simple);
%         flux_n(index) = Calc_RainBound(t,h(index),DELTAX(index),rain,prediction_data);
%         flux_n_old(index) = Calc_RainBound(t_old,h_old(index),DELTAX(index),rain,prediction_data);
%         flux_e(index) = -k_face(k_new,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H(index+1) - H(index))/dx(index));
%         flux_e_old(index) = -k_face(k_old,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H_old(index+1) - H_old(index))/dx(index));
%         
%         F(index) = psi_new(index) - psi_old(index) + ...
%             theta * dt * (1/DELTAX(index) * (flux_e(index) - ...
%             flux_w) + 1/DELTAZ(index) * (flux_n(index) - ...
%             flux_n(index-Nx)) - Q(index)) + ...
%             (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index) - ...
%             flux_w_old) + 1/DELTAZ(index) * (flux_n_old(index) - ...
%             flux_n_old(index-Nx)) - Q_old(index));
%     end
% end
% 
% for k = 1:num_indexes
%     index = index_list(k);
%     if i_list(k) == Nz && j_list(k) >=2 && j_list(k) <= Nx-1 % Top
%         flux_n(index) = Calc_RainBound(t,h(index),DELTAX(index),rain,prediction_data);
%         flux_n_old(index) = Calc_RainBound(t_old,h_old(index),DELTAX(index),rain,prediction_data);
%         flux_e(index) = -k_face(k_new,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H(index+1) - H(index))/dx(index));
%         flux_e_old(index) = -k_face(k_old,dx,index,1)* Kxz_face(Kxx,1,index,hetgen.boundary(index)) * ((H_old(index+1) - H_old(index))/dx(index));
%         F(index) = psi_new(index) - psi_old(index) + ...
%             theta * dt * (1/DELTAX(index) * (flux_e(index) - ...
%             flux_e(index-1)) + 1/DELTAZ(index) * (flux_n(index) - ...
%             flux_n(index-Nx)) - Q(index)) + ...
%             (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index) - ...
%             flux_e_old(index-1)) + 1/DELTAZ(index) * (flux_n_old(index) - ...
%             flux_n_old(index-Nx)) - Q_old(index));
%     end
% end
% 
% for k = 1:num_indexes
%     index = index_list(k);
%     if i_list(k) == Nz && j_list(k) == Nx % Top Right
%         flux_n(index) = Calc_RainBound(t,h(index),DELTAX(index),rain,prediction_data);
%         flux_n_old(index) = Calc_RainBound(t_old,h_old(index),DELTAX(index),rain,prediction_data);
%         
%         F(index) = psi_new(index) - psi_old(index) + ...
%             theta * dt * (1/DELTAX(index) * ( - ...
%             flux_e(index-1)) + 1/DELTAZ(index) * (flux_n(index) - ...
%             flux_n(index-Nx)) - Q(index)) + ...
%             (1 - theta) * dt * (1/DELTAX(index) * ( - ...
%             flux_e_old(index-1)) + 1/DELTAZ(index) * (flux_n_old(index) - ...
%             flux_n_old(index-Nx)) - Q_old(index));
%     end
% end
% 
% for k = 1:num_indexes
%     index = index_list(k);
%     if i_list(k)>=2 && i_list(k)<= Nz-1 && j_list(k) >=2 && j_list(k)<= Nx-1 % Internal
%         F(index) = psi_new(index) - psi_old(index) + ...
%             theta * dt * (1/DELTAX(index) * (flux_e(index) - ...
%             flux_e(index-1)) + 1/DELTAZ(index) * (flux_n(index) - ...
%             flux_n(index-Nx)) - Q(index)) + ...
%             (1 - theta) * dt * (1/DELTAX(index) * (flux_e_old(index) - ...
%             flux_e_old(index-1)) + 1/DELTAZ(index) * (flux_n_old(index) - ...
%             flux_n_old(index-Nx)) - Q_old(index));
%     end
% end
% F = F(index_list(num_indexes));


end