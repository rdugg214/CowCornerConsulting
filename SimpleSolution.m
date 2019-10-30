function [h] = SimpleSolution(SET)
v2struct(SET);
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
zonebin = zeros(length(x),4);
for i = 1:Nz
    for j = 1:Nx
        index = (i-1)*Nx + j;
        if ((0 <= x(index) && x(index) <= 50) && (40 <= z(index) && z(index) <= 100))...
                || ((50 < x(index) && x(index) <= 300) && (50 < z(index) && z(index) <= 100))
            zonetype(index) = 1;
            zonebin(index,1) = 1;
        end
        if  ((0 <= x(index) && x(index) <= 500) && (0 <= z(index) && z(index) <= 40))...
                || ((350 < x(index) && x(index) <= 500) && (40 < z(index) && z(index) <= 50))
            zonetype(index) = 3;
            zonebin(index,2) = 1;
        end
        if  ((50 < x(index) && x(index) <= 350) && (40 < z(index) && z(index) <= 50))
            zonetype(index) = 2;
            zonebin(index,3) = 1;
        end
        if  ((300 < x(index) && x(index) <= 500) && (50 <= z(index) && z(index) <= 100))
            zonetype(index)=4;
            zonebin(index,4) = 1;
        end
        end
end
zonebin = logical(zonebin);
hetgen.boundary = sum(zonebin,2)>1;
hetgen.Kxx = [2.6 0.08 3.9 4.3];
hetgen.Kzz = [0.91 0.016 1.17 0.21];
hetgen.psi_res = [0.01 0.106 0.01 0.01];
hetgen.psi_sat = [0.33 0.4686 0.1 0.075];
hetgen.alpha = [1.43 1.04 2.8 2.5];
hetgen.n = [1.51 1.3954 2.239 2.0];
hetgen.m = 1 - (1./hetgen.n);
hetgen.xcos = [1 -1 1 -1]; %NE NW SW SE
hetgen.zcos = [Nx Nx -Nx -Nx];
hetgen.xnos = [0 -1 -1 0]; %E N W S
hetgen.znos = [0 0 -Nx -Nx];
%% Sub in all.
alpha = zeros(size(x)); n = zeros(size(x)); m = zeros(size(x));
psi_res = zeros(size(x)); psi_sat = zeros(size(x));Kxx = zeros(size(x));
Kzz = zeros(size(x));
for i = 1:Nz
    for j = 1:Nx
        index = (i-1)*Nx + j;
         if hetgen.boundary(index)
             Kxx(index) = mean(hetgen.Kxx(zonebin(index,:))); 
            Kzz(index) = mean(hetgen.Kzz(zonebin(index,:)));
            psi_res(index) =  mean(hetgen.psi_res(zonebin(index,:)));
            psi_sat(index) = mean(hetgen.psi_sat(zonebin(index,:))); 
            alpha(index) = mean(hetgen.alpha(zonebin(index,:))); 
            n(index) =  mean(hetgen.n(zonebin(index,:)));
        else
            Kxx(index) = hetgen.Kxx(zonetype(index)); 
            Kzz(index) = hetgen.Kzz(zonetype(index));
            psi_res(index) = hetgen.psi_res(zonetype(index));
            psi_sat(index) =hetgen.psi_sat(zonetype(index)); 
            alpha(index) =hetgen.alpha(zonetype(index)); 
            n(index) = hetgen.n(zonetype(index));
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
% zind = 3;
% n(:) = n(zind);
% m(:) = m(zind);
% psi_res(:) = psi_res(zind);
% psi_sat(:) = psi_sat(zind);
% Kxx(:) = Kxx(zind);
% Kzz(:) = Kzz(zind);
% alpha(:) = alpha(zind);
load prediction_data.mat
params = {N, Nx, Nz, alpha, n , m, psi_res, psi_sat, x , z, ...
    Kxx , Kzz , dx , dz, DELTAX, DELTAZ,t_on_CSG,t_on_PUMP, simple,Pr,hetgen,prediction_data};

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
hload = load('h_init.mat');
h_old = hload.h;

%% LOOP!
% Define static variables
t =0;
t_old = 0;
h = h_old;
figm = figure('units','normalized','outerposition',[0.2 0.3 0.8 0.7]);
omega = 1;
dt = dt_o;
t_hist = 0;

%---- plot average vs current
psi_now = helper_getpsinow(h_old, alpha,n,m,psi_res,psi_sat,x,z,dx,dz,hetgen);
psi_av  = 1/(Lx*Lz) * sum(DELTAX.*DELTAZ.*psi_now);
psi_av_hist = psi_av;
psi_int = psi_av;
Tarf =0.6544;
psi_guess_func = @(t) psi_int + -(t.*(-1.803e-3)-sin(t.*pi.*(2.0./3.65e+2)).*1.047390722740609e-1)/Lz...
    - ( ((Tarf * Pr)/((365)*(75-55)) ) *t )/(Lx*Lz);
psi_guess =  psi_guess_func(t);
psi_guess_hist = psi_av;
%----- capture outputs from flux
% riverloc = x == 0 & z>=80 & z<=100;
% CSGloc = x == 500 & z <= 5 & z >=0;
 riverloc = z >=80 & x ==0 ;
 CSGloc = z <=5 & x ==500;
rainloc = z ==100;
Ballocs = [riverloc';CSGloc';rainloc'];
Balvals = zeros(size(Ballocs,1),1);
Ballarr = Balvals;
% Ballarr = 0;

helper_plotcmap(X,Z,helper_row2mat(Nz,Nx,zonetype),helper_row2mat(Nz,Nx,zonetype),figm);

if SAVEVID
movegui(figm,'onscreen');
dname = char(datetime);
dname(dname==' ') = '_'; dname = ['_',dname,'_'];
dname(dname == ':') = '-';
vidObj = VideoWriter([dname '.avi']);
vidObj.Quality = 100;
vidObj.FrameRate = 1;
open(vidObj);
end
while t<endtime
    ftsuccess = true;
    success = false;
   
    while success == false
     dt = omega*dt;
    t = t_old+dt;
    F = @(h) FVM(h,dt, t, t_old,h_old, params);
    Jacobian = @(F,x,Fx0) NEW_JacobianFD(F,x,Fx0,params,dt, t, t_old, h_old);
    [h,success] = NEW_Newton_Solver(F,h_old,Jacobian, "Shamanskii");
    [~,hgain] =  FVM(h,dt, t, t_old,h_old, params);
    if success == false
     
        ftsuccess = false;
    fprintf('New dt = %3.2f\n',dt);
    end
    omega =0.5;
    end
    if ftsuccess == true && (dt/omega)/omega <dtmax
    dt = min(dtmax,(dt/omega)/omega);
    fprintf('New dt = %3.8f\n',dt);
    end
%     if ftsuccess == true & (dt/omega)/omega <dtmax
%     dt = min(dtmax,(dt/omega)/omega);
%      fprintf('New dt = %3.8f\n',dt);
%     end
    fprintf('Success!\n')
    
    t_hist = [t_hist t];
    %
%     Balvals = Balvals + [sum(hgain.w(Ballocs(1,:)).*DELTAZ(Ballocs(1,:))); ...
%     sum(hgain.e(Ballocs(2,:)).*DELTAZ(Ballocs(2,:))) ];
    Balvals = [sum(hgain.w' ./DELTAX(Ballocs(1,:))); ...
    sum(hgain.e(Ballocs(2,:))./DELTAX(Ballocs(2,:)));...
     sum(hgain.n(Ballocs(3,:)) ./DELTAZ(Ballocs(3,:)))];
%     Ballarr(:,end+1) = Ballarr(:,end) +  Balvals;
     Ballarr(:,end+1) =   Balvals;
    % comparison average vs current
    psi_old = psi_now;
    psi_now = helper_getpsinow(h, alpha,n,m,psi_res,psi_sat,x,z,dx,dz,hetgen);
    psi_av  = 1/(Lx*Lz) * sum(DELTAX.*DELTAZ.*psi_now);
    psi_av_hist = [psi_av_hist psi_av];
    psi_guess =   psi_guess_func(t);
    psi_guess_hist = [psi_guess_hist psi_guess];
   helper_plot_h_psi_av(Nz,Nx,figm,X,Z,psi_now,h,t_hist,psi_av_hist,psi_guess_hist,Ballarr)
   
   title(sprintf('t = %.2f (years) (current dt = %.2f (days))',t/365,dt));
   h_old = h;
   t_old = t;
   t = t+dt;
   drawnow
    if SAVEVID
        writeVideo(vidObj,getframe(gcf));
    end
end
if SAVEVID
 close(vidObj);
end
end