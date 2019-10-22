function [] = SimpleSolution(dtmax, endtime, t_on_CSG,t_on_PUMP,geometric)
%% Simiple version of the problem
% Solve the advection diffusion problem
% Problem setup
Lx = 500;
Lz = 100;
dt_o = 2;
%% Generate a uniform mesh with N node points
if length(geometric) == 3
    % Complex Mesh with each regions modelled seperately. 

% x1 = GP_sym(0,50,50/2,1.01); x2 = GP_sym(50,300,250/2,1.02); x3 = GP_sym(300,350,50/2,1.02);  x4 = GP_sym(350,500,150/2,1.02);
% 
% z1 = GP_sym(0,5,10/2,1.05);  z2 = GP_sym(5,40,35/2,1.05); z3 = GP_sym(40,50,10/2,1.03); z4 = GP_sym(50,80,30/2,1.05); z5 = GP_sym(80,100,20/2,1.05); 

% x=[x1(1:end-1),x2(1:end-1),x3(1:end-1),x4]; z = [z1(1:end-1),z2(1:end-1),z3(1:end-1),z4(1:end-1),z5];

% Simple Mesh with geometric progression
nx = geometric(2);
nz = geometri(3);
x = round(GP_sym(0,500,nx,1.3),4); z = round(GP_sym(0,100,ny,1.5),4);
o_x = x; o_z = z;
else 
    dx = geometric(1);
    dz = geometric(2);
o_x = 0:dx:Lx;
o_z = 0:dz:Lz;
end
Nx = length(o_x);
Nz = length(o_z);

N = Nx*Nz;
x = (zeros(N,1));
z = x;

dx = zeros(size(x));
dz = zeros(size(z));
DELTAX = zeros(size(x));
DELTAZ= zeros(size(z));
 for i = 1:Nz %loop by row (left to right, dowtn to top)
    for j = 1:Nx
        index = (i-1)*Nx + j;
        x(index) = o_x(j);
        z(index) = o_z(i);
        node_num_MAT(i,j) = index;
       node_num(index) = index;
    end
end

X = helper_row2mat(Nz,Nx,x);
Z =helper_row2mat(Nz,Nx,z);

for i = 1:Nz %calc dxdz
    for j = 1:Nx
        index = (i-1)*Nx + j;
        if j == 1
            dx(index) =(x(index+1)-x(index));
        elseif j == Nx
            dx(index) =(x(index)-x(index-1));
        else
            dx(index) = x(index) - x(index-1);
        end
        
        if i == 1
            dz(index) = (z(index+Nx)-z(index));
        elseif i == Nz
            dz(index) = (z(index)-z(index-Nx));
        else
            dz(index) =z(index) -  z(index-Nx);
        end  
    end
end

for i = 1:Nz
    for j = 1:Nx
        index = (i-1)*Nx + j;
        if j == 1
            DELTAX(index) = (x(index+1)-x(index))/2;
        elseif j == Nx
            DELTAX(index) = (x(index)-x(index-1))/2;
        else
            DELTAX(index) = (dx(index-1)+dx(index))/2;
        end
        if i == 1
            DELTAZ(index) = (z(index+Nx)-z(index))/2;
        elseif i == Nz
            DELTAZ(index) = (z(index)-z(index-Nx))/2;
        else
            DELTAZ(index) = (dz(index-1)+dz(index))/2;
        end
    end
end %Calc DXDZ
%% CALCULATE ZONE TYPES

% Generate params %NO OVERLAP FOR NOW
zonetype = zeros(size(x));
for i = 1:Nz
    for j = 1:Nx
        index = (i-1)*Nx + j;
        if ((0 <= x(index) && x(index) <= 50) && (40 <= z(index) && z(index) <= 100))...
                || ((50 < x(index) && x(index) <= 300) && (50 < z(index) && z(index) <= 100))
            zonetype(index) = 1;
        elseif  ((0 <= x(index) && x(index) <= 500) && (0 <= z(index) && z(index) <= 40))...
                || ((350 < x(index) && x(index) <= 500) && (40 < z(index) && z(index) <= 50))
            zonetype(index) = 2;
        elseif  ((50 < x(index) && x(index) <= 350) && (40 < z(index) && z(index) <= 50))
            zonetype(index) = 3;
        elseif  ((300 < x(index) && x(index) <= 500) && (50 <= z(index) && z(index) <= 100))
            zonetype(index)=4;
        end
        
    end
end



%% Sub in all.
alpha = zeros(size(x)); n = zeros(size(x)); m = zeros(size(x));
psi_res = zeros(size(x)); psi_sat = zeros(size(x));Kxx = zeros(size(x));
Kzz = zeros(size(x));
for i = 1:Nz
    for j = 1:Nx
        index = (i-1)*Nx + j;
        switch zonetype(index)
            case 1
                Kxx(index) =2.6; Kzz(index) = 0.91; psi_res(index) = 0.01;
                psi_sat(index) = 0.33; alpha(index) = 1.43; n(index) = 1.51;
            case 3
                Kxx(index) =0.08; Kzz(index) = 0.016; psi_res(index) = 0.106;
                psi_sat(index) = 0.4686; alpha(index) = 1.04; n(index) = 1.3954;
            case 2
                Kxx(index) =3.9; Kzz(index) = 1.17; psi_res(index) = 0.01;
                psi_sat(index) = 0.1; alpha(index) = 2.8; n(index) = 2.239;
            case 4
                Kxx(index) =4.3; Kzz(index) = 0.21; psi_res(index) = 0.01;
                psi_sat(index) = 0.075; alpha(index) = 2.5; n(index) = 2.0;
        end
        
    end
end
m = 1 - (1./n);

n(:) = mean(n);
m(:) = mean(m);
psi_res(:) = mean(psi_res);
psi_sat(:) = mean(psi_sat);
Kxx(:) = mean(Kxx);
Kzz(:) = mean(Kzz);
alpha(:) = mean(alpha);
% n(:) = n(2);
% m(:) = m(2);
% psi_res(:) = psi_res(2);
% psi_sat(:) = psi_sat(2);
% Kxx(:) = Kxx(2);
% Kzz(:) = Kzz(2);
% alpha(:) = alpha(2);
params = {N, Nx, Nz, alpha, n , m, psi_res, psi_sat, x , z, ...
    Kxx , Kzz , dx , dz, DELTAX, DELTAZ,t_on_CSG,t_on_PUMP};

%% Generate Initial Solution
h_old = zeros(size(x));
hbot = -5;
htop = -10;
for i = 1:Nz
    for j = 1:Nx
        index = (i-1)*Nx + j;
        h_old(index) = hbot + ((htop - hbot)*z(index)/Lz);
    end
end

%% LOOP!
% Define static variables
t =0;
t_old = 0;
figure;
h = h_old;
figm = figure('units','normalized','outerposition',[0.2 0.3 0.8 0.7]);
omega = 1;
dt = dt_o;
t_hist = 0;

psi_now = helper_getpsinow(h_old, alpha,n,m,psi_res,psi_sat);
psi_av  = 1/(Lx*Lz) * sum(DELTAX.*DELTAZ.*psi_now);
psi_av_hist = psi_av;
psi_int = psi_av;
Tarf =0.6544;
Pr = 0.25;
psi_guess_func = @(t) psi_int + -(t.*(-1.803e-3)-sin(t.*pi.*(2.0./3.65e+2)).*1.047390722740609e-1)/Lz...
    - ( ((Tarf * Pr)/((365)*(75-55)) ) *t )/(Lx*Lz);
psi_guess =  psi_guess_func(t);
psi_guess_hist = psi_av;
helper_plot_h_psi(Nz,Nx,figm,X,Z,psi_now,h)


    drawnow
while t<endtime
    ftsuccess = true;
    success = false;
    tic
   
    while success == false
    dt = omega*dt;
    t = t_old+dt;
%     h = fsolve(@(h) FVM(h,dt, t, t_old,h_old, params),h_old,optimoptions('fsolve','Display','off')); 
%     tic
    F = @(h) FVM(h,dt, t, t_old,h_old, params);
    F_in = @(h, index) FVM_index(h, dt, t, t_old, h_old, params, index);
%     toc
    Jacobian = @(F,x,Fx0,F_in) NEW_JacobianFD(F,x,Fx0,F_in);
    [h,success] = NEW_Newton_Solver(F,h_old,Jacobian, F_in);
    if success == false;
        ftsuccess = false;
    fprintf('New dt = %3.2f\n',dt);
    end
    omega =0.5;
    end
    if ftsuccess == true & (dt/omega)/omega <dtmax
    dt = min(dtmax,(dt/omega)/omega);
     fprintf('New dt = %3.8f\n',dt);
    end
    fprintf('Success!\n')
%     figure 
    psi_old = psi_now;
    psi_now = helper_getpsinow(h, alpha,n,m,psi_res,psi_sat);
    h_old = h;
    t_old = t;
    
    psi_av  = 1/(Lx*Lz) * sum(DELTAX.*DELTAZ.*psi_now);
    t_hist = [t_hist t];
    psi_av_hist = [psi_av_hist psi_av];
    psi_guess =   psi_guess_func(t);
    psi_guess_hist = [psi_guess_hist psi_guess];
    
   
   helper_plot_h_psi_av(Nz,Nx,figm,X,Z,psi_now./psi_sat,h,t_hist,psi_av_hist,psi_guess_hist,1)
   title(sprintf('t = %.2f (years) (current dt = %.2f (days))',t/365,dt));
   t = t+dt;
    drawnow
end
end