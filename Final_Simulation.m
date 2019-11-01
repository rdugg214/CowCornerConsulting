function [RES] = Final_Solution(SET)
v2struct(SET);
%% Produce a solution and simulation of the problem

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIXED PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Environment setup
Lx = 500; %Length in x
Lz = 100; %Lenght in z
dt_o = 3; %Original timestep

%% Generate a Mesh
if length(geometric) == 3 %If length of geometric properties given is 3, then apply geometric progressoin
    
    % Complex Mesh with each regions modelled seperately.
    nx = geometric(2); %number of (desired) nodes in x
    nz = geometric(3); %number of (desired) nodes in z
    
    % Create geometric progression in x components for each region, scaled by
    % size of region
    x1 = GP_sym(0,50,ceil(nx/10),1.01); x2 = GP_sym(50,300,ceil(5*nx/10),1.02); x3 = GP_sym(300,350,ceil(3*nx/10),1.02);  x4 = GP_sym(350,500,ceil(1*nx/10),1.02);
    % Create geometric progression in z components for each region
    z1 = GP_sym(0,5,ceil(nz/10),1.05);  z2 = GP_sym(5,40,ceil(3*nz/10),1.05); z3 = GP_sym(40,50,ceil(nz/10),1.03); z4 = GP_sym(50,80,ceil(3*nz/10),1.05); z5 = GP_sym(80,100,ceil(2*nz/10),1.05);
    
    %concatenate all x values
    x=[x1(1:end-1),x2(1:end-1),x3(1:end-1),x4]; z = [z1(1:end-1),z2(1:end-1),z3(1:end-1),z4(1:end-1),z5];
    
    % set this as the 'reference xz values'
    o_x = x; o_z = z;
    
else %otherwise, if length is 2, given is fixed dz and create uniform
    dx = geometric(1);
    dz = geometric(2);
    
    % create base x and z values
    o_x = 0:dx:Lx; o_z = 0:dz:Lz;
end

% Gather the number of nodes for each axis
Nx = length(o_x); Nz = length(o_z);
% Then gather the number of nodes together
N = Nx*Nz;

% Initialize the vector of x,dx, Delta x and z,dz, Delta z for each node
x = (zeros(N,1)); z = x;
dx = zeros(size(x)); dz = zeros(size(z));
DELTAX = zeros(size(x)); DELTAZ= zeros(size(z));

% Enter x and z values
for i = 1:Nz %loop by row (left to right, bottom to top)
    for j = 1:Nx
        index = (i-1)*Nx + j;
        x(index) = o_x(j);
        z(index) = o_z(i);
        node_num_MAT(i,j) = index;
        node_num(index) = index;
    end
end

% Create matrix version of x and z row vectors for mesh
X = helper_row2mat(Nz,Nx,x);
Z =helper_row2mat(Nz,Nx,z);

% Enter dx and dz values
for i = 1:Nz
    for j = 1:Nx
        index = (i-1)*Nx + j;
        if j == 1 %if left boundary
            dx(index) =(x(index+1)-x(index));
        elseif j == Nx %if right boundary
            dx(index) =(x(index)-x(index-1));
        else %if middle
            dx(index) = x(index) - x(index-1); %calculate spacing between node
        end
        
        if i == 1 %if bottom boundary
            dz(index) = (z(index+Nx)-z(index));
        elseif i == Nz %if top boundary
            dz(index) = (z(index)-z(index-Nx));
        else %if middle
            dz(index) =z(index) -  z(index-Nx); %calculate spacing between node
        end
    end
end
% Enter Dx and Dz values
for i = 1:Nz
    for j = 1:Nx
        index = (i-1)*Nx + j;
        if j == 1 %if left boundary
            DELTAX(index) = (x(index+1)-x(index))/2;
        elseif j == Nx %if right boundary
            DELTAX(index) = (x(index)-x(index-1))/2; %calculate width of  CV (1/2 of middle)
        else  %if middle
            DELTAX(index) = (dx(index-1)+dx(index))/2;  %calculate width of CV
        end
        if i == 1 %if bottom boundary
            DELTAZ(index) = (z(index+Nx)-z(index))/2;
        elseif i == Nz %if top boundary
            DELTAZ(index) = (z(index)-z(index-Nx))/2;  %calculate height of  CV (1/2 of middle)
        else  %if middle
            DELTAZ(index) = (dz(index-1)+dz(index))/2;  %calculate height of CV
        end
    end
end %Calc DXDZ
%% Allocate Zone Property Information
% Singular value zone type for each node (ith entry = kth zone), for
% homogenous case
zonetype = zeros(size(x));
% Logical value zone type (size N x K) for each node (ith entry = 1 or 0
% for kth zone) for heterogenous case
zonebin = zeros(length(x),4);

for i = 1:Nz
    for j = 1:Nx
        index = (i-1)*Nx + j;
        % Alluvium
        if ((0 <= x(index) && x(index) <= 50) && (40 <= z(index) && z(index) <= 100))...
                || ((50 < x(index) && x(index) <= 300) && (50 < z(index) && z(index) <= 100))
            zonetype(index) = 1;
            zonebin(index,1) = 1;
        end
        % Walloon Coal
        if  ((0 <= x(index) && x(index) <= 500) && (0 <= z(index) && z(index) <= 40))...
                || ((350 < x(index) && x(index) <= 500) && (40 < z(index) && z(index) <= 50))
            zonetype(index) = 3;
            zonebin(index,3) = 1;
        end
        % Confining Layer
        if  ((50 < x(index) && x(index) <= 350) && (40 < z(index) && z(index) <= 50))
            zonetype(index) = 2;
            zonebin(index,2) = 1;
        end
        % Main range volcanics
        if  ((300 < x(index) && x(index) <= 500) && (50 <= z(index) && z(index) <= 100))
            zonetype(index)=4;
            zonebin(index,4) = 1;
        end
    end
end

zonebin = logical(zonebin); %convert to logical

% Heterogenous related variables
hetgen.boundary = sum(zonebin,2)>1; %for any nodes with multiple properties detected, label as boundary
hetgen.Kxx = [2.6 0.08 3.9 4.3]; %store the values according to k index
hetgen.Kzz = [0.91 0.016 1.17 0.21];
hetgen.psi_res = [0.01 0.106 0.01 0.01];
hetgen.psi_sat = [0.33 0.4686 0.1 0.075];
hetgen.alpha = [1.43 1.04 2.8 2.5];
hetgen.n = [1.51 1.3954 2.239 2.0];
hetgen.m = 1 - (1./hetgen.n);

% for heterogeneous, attach the offsets to access the corners (xcos,zcos)
% and the nodes (xnos,znos) directly.
hetgen.xcos = [1 -1 1 -1]; %NE NW SW SE
hetgen.zcos = [Nx Nx -Nx -Nx];
hetgen.xnos = [0 -1 -1 0]; %E N W S
hetgen.znos = [0 0 -Nx -Nx];

% For example to access the NORTHEAST node from the current index then we
% use: x(index+hetgen.xcos(1)+hetgen.zcos(1))
% To access the SOUTHWEST node then:
% use: x(index+hetgen.xcos(3)+hetgen.zcos(3))
% To access the NORTH node then:
% use: x(index+hetgen.xnos(2)+hetgen.znos(2))
% To access the S node then:
% use: x(index+hetgen.xnos(4)+hetgen.znos(4))

% Create the properties vectors for each node:
alpha = zeros(size(x)); n = zeros(size(x)); m = zeros(size(x));
psi_res = zeros(size(x)); psi_sat = zeros(size(x));Kxx = zeros(size(x));
Kzz = zeros(size(x));
for i = 1:Nz
    for j = 1:Nx
        index = (i-1)*Nx + j;
        %          if hetgen.boundary(index) %Heterogenous case, if node is a
        %          %boundary, take the mean of the values
        %              Kxx(index) = mean(hetgen.Kxx(zonebin(index,:)));
        %             Kzz(index) = mean(hetgen.Kzz(zonebin(index,:)));
        %             psi_res(index) =  mean(hetgen.psi_res(zonebin(index,:)));
        %             psi_sat(index) = mean(hetgen.psi_sat(zonebin(index,:)));
        %             alpha(index) = mean(hetgen.alpha(zonebin(index,:)));
        %             n(index) =  mean(hetgen.n(zonebin(index,:)));
        %         else
        % For other nodes, just take the singular type they were
        % assigned in
        Kxx(index) = hetgen.Kxx(zonetype(index));
        Kzz(index) = hetgen.Kzz(zonetype(index));
        psi_res(index) = hetgen.psi_res(zonetype(index));
        psi_sat(index) =hetgen.psi_sat(zonetype(index));
        alpha(index) =hetgen.alpha(zonetype(index));
        n(index) = hetgen.n(zonetype(index));
        
        %For semi-heterogenous case, we do not take the mean since this
        %produces incorrect values when not entirety of heterogenous is
        %implemented
        %         end
        
    end
end
m = 1 - (1./n);

% % Truly homogenous case (same properties THROUGHOUT the aquifer), then
% % assign the properties for all nodes as one, mean value
% n(:) = mean(n);
% m(:) = mean(m);
% psi_res(:) = mean(psi_res);
% psi_sat(:) = mean(psi_sat);
% Kxx(:) = mean(Kxx);
% Kzz(:) = mean(Kzz);
% alpha(:) = mean(alpha);

% if drought type rain, assign the drought data, which is the regular
% prediction data, but with any entry <620mm/year reduced to 5% of the
% predicted, based on severe drought classification
% http://www.bom.gov.au/climate/glossary/drought.shtml
load('prediction_data.mat')
if rain == 3
    prediction_data(prediction_data<620) = prediction_data(prediction_data<620)*0.05;
end


params = {N, Nx, Nz, alpha, n , m, psi_res, psi_sat, x , z, ...
    Kxx , Kzz , dx , dz, DELTAX, DELTAZ,t_on_CSG,t_on_PUMP, simple,Pr,hetgen,prediction_data,DELCSG,rain};

%% Generate Initial Solution
h_old = zeros(size(x));
if length(h_init) == 0 %Standard initial values
    hbot = -5;
    htop = -10;
    
    for i = 1:Nz
        for j = 1:Nx
            index = (i-1)*Nx + j;
            h_old(index) = hbot + ((htop - hbot)*z(index)/Lz);
        end
    end
elseif length(h_init) == 1 %Saturated initial values
    hbot = 95;
    htop = -2;
    for i = 1:Nz
        for j = 1:Nx
            index = (i-1)*Nx + j;
            h_old(index) = hbot + ((htop - hbot)*z(index)/Lz);
        end
    end
else
    h_old = h_init;
end

%% LOOP!

% Initialise t and t_old
t =t_init; t_old = t_init;
h = h_old;

% Initialise figures for water profile (figm) and fluxes (figb)
figm = figure('units','normalized','outerposition',[0.2 0.5 0.8 0.5]);
figb = figure('units','normalized','outerposition',[0.2 0.05 0.8 0.4]);

% initialise omega parameter for varying step size
omega = 1;
dt = dt_o; %initial step size

%---- plot average vs current
% get the psi vector for current time
psi_now = helper_getpsinow(h_old, alpha,n,m,psi_res,psi_sat,x,z,dx,dz,hetgen);
% get average psi value for system
psi_av  = 1/(Lx*Lz) * sum(DELTAX.*DELTAZ.*psi_now);
%store as a vector across time
psi_av_hist = psi_av;

psi_int = psi_av; %initial psi value for `analytical guess'
Tarf =0.6544; %Average rainfall for simple rain
%Guess function including rainfall and pumpin
psi_guess_func = @(t) psi_int + -(t.*(-1.803e-3)-sin(t.*pi.*(2.0./3.65e+2)).*1.047390722740609e-1)/Lz...
    - ( ((Tarf * Pr)/((365)*(75-55)) ) *t )/(Lx*Lz);
%gather guess at current time
psi_guess =  psi_guess_func(t);
%store as vector
psi_guess_hist = psi_av;

%----- capture outputs from flux
% Gather locations of desired outputs
riverloc = z >=80 & x ==0 ;
CSGloc = z <=5 & x ==500;
rainloc = z ==100;
evaploc = ((50<=x & x <= 100) & (85 <= z & z<= 100))...
    |  ((100<=x & x <= 300) & (95 <= z & z<= 100)) ...
    | ((300<=x & x <= 500) & (90 <= z & z<= 100));
pumploc = x==100  & 55 <= z & z<= 75;

% store as a matrix
Ballocs = [riverloc';CSGloc';rainloc';evaploc';pumploc'];
% initialize values matrix
Balvals = zeros(size(Ballocs,1),1);
Ballarr = Balvals;

%gather average total pressure head
Harr =1/(Lx*Lz) * sum(h_old+z);

% initialise time vector
t_vector = [t];
H_vector = [mean(h+z)];
if SAVEVID
    % Store video of outputs
    movegui(figb,'onscreen');
    dname = ['_',NAME,'_','GRAPH'];
    vidObj = VideoWriter([dname '.avi']);
    vidObj.Quality = 100;
    vidObj.FrameRate = 1;
    open(vidObj);
    % Store video of water profile
    movegui(figm,'onscreen');
    dname = ['_',NAME,'_','CMAP'];
    vidObj2 = VideoWriter([dname '.avi']);
    vidObj2.Quality = 100;
    vidObj2.FrameRate = 1;
    open(vidObj2);
end

%While simulation time is met
while t<endtime
    ftsuccess = true; %initialise test for `first time success'
    success = false; %initialise test for `success' of solver
    
    while success == false
        % update delta t
        dt = omega*dt;
        t = t_old+dt;
        [k_old, psi_old, Q_old, ~] = FVM_pre_calcs(h_old, dt, t, params);
        F = @(h) FVM(h,dt, t, t_old,h_old, k_old, psi_old, Q_old, params);
        Jacobian = @(F,x,Fx0) NEW_JacobianFD(F,x,Fx0,params,dt, t, t_old, h_old);
        [h,success] = NEW_Newton_Solver(F,h_old,Jacobian, "Shamanskii");
        [~,hgain] =  FVM(h,dt, t, t_old,h_old, k_old, psi_old, Q_old, params);
        if success == false
            ftsuccess = false;
            fprintf('Current dt = %3.2f\n',dt);
        end
        omega =0.5;
    end
    H_vector(length(H_vector)+1) = mean(h+z);
    t_vector(length(t_vector)+1) = t;
    if ftsuccess == true && (dt/omega)/omega <dtmax
        dt = min(dtmax,(dt/omega)/omega);
        fprintf('Current dt = %3.8f\n',dt);
    else
        fprintf('Current dt = %3.8f\n',dt);
    end
    %     if ftsuccess == true & (dt/omega)/omega <dtmax
    %     dt = min(dtmax,(dt/omega)/omega);
    %      fprintf('New dt = %3.8f\n',dt);
    %     end
    fprintf('Success!\n')
    

    Balvals = [sum(hgain.w'); ...
        sum(hgain.e(Ballocs(2,:)));...
        sum(hgain.n(Ballocs(3,:)));...
        sum(hgain.Q(Ballocs(4,:)));...
        sum(hgain.Q(Ballocs(5,:)))];
    %     Ballarr(:,end+1) = Ballarr(:,end) +  Balvals;
    Ballarr(:,end+1) =   Balvals;
    Harr(1,end+1) = 1/(Lx*Lz) * sum(h+z);
    % comparison average vs current
    psi_old = psi_now;
    psi_now = helper_getpsinow(h, alpha,n,m,psi_res,psi_sat,x,z,dx,dz,hetgen);
    psi_av  = 1/(Lx*Lz) * sum(DELTAX.*DELTAZ.*psi_now);
    psi_av_hist = [psi_av_hist psi_av];
    psi_guess =   psi_guess_func(t);
    psi_guess_hist = [psi_guess_hist psi_guess];
    helper_plot_h_psi_av(Nz,Nx,figm,X,Z,psi_now./psi_sat,h,t_vector,psi_av_hist,psi_guess_hist,t,dt,simple,DELCSG,Pr)
    helper_plot_Ballarr(t_vector,Ballarr,figb,DELCSG,Pr,t,dt);
    
    h_old = h;
    t_old = t;
    t = t+dt;
    drawnow
    if SAVEVID
        writeVideo(vidObj,getframe(figb));
        writeVideo(vidObj2,getframe(figm));
    end
end
% figure;
% method_times = [method_times(1, :), method_times(2, :)];
% b = bar(method_times);
% b.FaceColor = 'flat';
% b.CData(1:3,:) = [1 0 0; 1 0 0; 1 0 0];
% set(gca,'xticklabel',{"Blackslash Full Newton", "Blackslash Chord", "Blackslash Shamanskii", "GMRES Full Newton", "GMRES Chord", "GMRES Shamanskii"});
% title("Time to execute by newton solver method");
% ylabel("Time (s)");
% xtickangle(45);
if SAVEVID
    close(vidObj);
    close(vidObj2);
end
RES.h_final = h;
RES.Ballarr = Ballarr;
RES.psi_av_hist = psi_av_hist;
RES.t_final = t;
RES.t_on_CSG = t_on_CSG;
RES.t_on_PUMP = t_on_PUMP;
RES.DELCSG = DELCSG;
RES.Pr =Pr;
RES.Harr = Harr;
RES.geometric = geometric;
RES.t_vector = t_vector;
RES.H_vector = H_vector;

end