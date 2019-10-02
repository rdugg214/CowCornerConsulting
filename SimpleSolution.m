%% Simiple version of the problem
% Solve the advection diffusion problem
clear; clc; close all;
% Problem setup
Lx = 500;
Lz = 100;


%% Generate a uniform mesh with N node points
dx = 10; dz = 10;

x = 0:dx:Lx;
Nx = length(x);
z = [];
Nz = length(0:dz:Lz);
N = Nx*Nz;
for i = 1:Nz
    z = [z, repmat((i-1)*dx,1,Nx)];
end


x = repmat(x, 1, Nz);
dx = zeros(size(x));
dz = zeros(size(z));
 [X,Z] = meshgrid(linspace(0,Lx,Nx),linspace(0,Lz,Nz));

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


% Generate params %NO OVERLAP FOR NOW
zonetype = zeros(size(x));
I = zeros(Nz,Nx,3);
for i = 1:Nz
    for j = 1:Nx
        index = (i-1)*Nx + j;
        if ((0 <= x(index) && x(index) <= 50) && (40 <= z(index) && z(index) <= 100))...
                || ((50 <= x(index) && x(index) <= 300) && (50 <= z(index) && z(index) <= 100))
            zonetype(index) = 1;
        elseif  ((0 <= x(index) && x(index) <= 500) && (0 <= z(index) && z(index) <= 40))...
                || ((350 <= x(index) && x(index) <= 500) && (40 <= z(index) && z(index) <= 50))
            zonetype(index) = 2;
        elseif  ((50 <= x(index) && x(index) <= 350) && (40 <= z(index) && z(index) <= 50))
            zonetype(index) = 3;
        elseif  ((300 <= x(index) && x(index) <= 500) && (50 <= z(index) && z(index) <= 100))
            zonetype(index)=4;
        end
        
    end
end

% for i = 1:Nz
%     for j = 1:Nx
%         index = (i-1)*Nx + j;
%         switch zonetype(index)
%             case 1
%                 I(i,j,1) = 255;
%             case 2
%                 I(i,j,2) = 255;
%             case 3
%                 I(i,j,3) = 255;
%             case 4
%                 I(i,j,1) = 255; I(i,j,2) = 255; I(i,j,3) = 255;
%         end
%     end
% end
% figure
% imshow(I)

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

 MAT_dx = helper_row2mat(Nz,Nx,h_old) ;
helper_plotcmap(X,Z,MAT_dx)
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
            case 2
                Kxx(index) =0.08; Kzz(index) = 0.016; psi_res(index) = 0.106;
                psi_sat(index) = 0.4686; alpha(index) = 1.04; n(index) = 1.3954;
            case 3
                Kxx(index) =3.9; Kzz(index) = 1.18; psi_res(index) = 0.01;
                psi_sat(index) = 0.1; alpha(index) = 2.8; n(index) = 2.239;
            case 4
                Kxx(index) =4.3; Kzz(index) = 0.21; psi_res(index) = 0.01;
                psi_sat(index) = 0.075; alpha(index) = 2.5; n(index) = 2.0;
        end
        
    end
end
m = 1 - (1./n);
params = {N, Nx, Nz, alpha, n , m, psi_res, psi_sat, x , z, ...
    Kxx , Kzz , dx , dz, DELTAX, DELTAZ,alpha};

%% LOOP!
% Define static variables
dt = 10; endtime = 100000;
figure;
t = 0+dt;
t_old =t;

% tic
% FVM(h_old,dt, t, t_old,h_old, params)
% toc

% h = NewtonSolver(h_old, 0);
% F = FVM(0,dt, t, t_old,h_old, params);
% fim = imag(F);
% [X,Y] = meshgrid(linspace(0,Lx,Nx),linspace(0,Lz,Nz));
% Z = zeros(size(Y));
% Y = flip(Y);
%  for i = 1:Nz
%     for j = 1:Nx
%         index = (i-1)*Nx + j;
%         Z(i,j) = fim(index);
%     end
% end
% surf(X,Y,Z)
% disp('err?')
while t<endtime
    h = fsolve(@(h) FVM(h,dt, t, t_old,h_old, params),h_old,optimoptions('fsolve','Display','off','FunctionTolerance',1e-1));
    t_old = t;
    t = t+dt;
    %     scatter(x,h);
 MAT_dx = helper_row2mat(Nz,Nx,h) ;
helper_plotcmap(X,Z,MAT_dx)
   title(['t = ' num2str(t)])
    h_old = h;
    drawnow
 
    %     h = NewtonSolver(h, t);
end

%% Newton Solver
% h_new = fsolve(@(h) FVM(h,dt, t, t_old,h_old, params),h_old,optimoptions('fsolve','Display','off','FunctionTolerance',1e-1));
%  MAT_dx = helper_row2mat(Nz,Nx,h_new) ;
% helper_plotcmap(X,Z,MAT_dx)
NewtonSolver(h_old, h_old, t, t_old, dt, params)