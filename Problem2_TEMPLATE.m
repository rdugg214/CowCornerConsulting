%% Benchmark Problem 1.2
% Solve the advection diffusion problem

clear; % Clear all variables from the current workspace
clc; % clear the command window
close all;
% It is important to clear before each test so variables from previous
% scripts/tests don't conflict with your current test

% Problem setup
U0 = 1; % Concentration at Boundary
D = 2.5e-3; % diffusion coefficient
V = 0.5; % Advective velocity
lambda = 0.05; % Source Term coefficient 
L = 1; % Length of domain
N = 51; % Number of node points on mesh
dt = 2^-3; % timestep size
endtime = 200; % End time
theta = 1; % Temporal weightings, 1 - Backward Euler, 0.5 - Crank-Nicholson, 0 - Forward Euler
sigma = 0; % Upwinding = 0, Downwind = 2, Centred = 1;
plot_times = 0:2^-2:endtime; % vector of times to save and plot the numerical solution
realtime_plot = true; % Should a plot of the solution be produced in realtime?

% You can clean up the number of inputs to a function by storing multiple
% parameters in a vector
parms = [theta L N D V lambda U0 sigma];

% Specify the boundary conditions
% Left Boundary: A0*U - B0*dU/dx = C0
A0 = 1; B0 = 0; C0 = U0;
% Right Boundary: AL*U + BL*dU/dx = CL
AL = 1; BL = 0; CL = 0;

% Put all Boundary condition parameters into an array
BC = [A0 B0 C0; AL BL CL];

% Generate a uniform mesh with N node points
x = linspace(0,L,N);

% Calculate Mesh Spacings
deltax = zeros(1,N-1); % distances between nodes
DELTAX = zeros(1,N); % distances between control volume faces
deltax(1) = x(2) - x(1); % Difference between first 2 nodes
deltax(N-1) = x(N)-x(N-1); % Difference between last 2 nodes
DELTAX(1) = (deltax(1))/2; % left boundary
DELTAX(N) = (deltax(N-1))/2; % right boundary
% internal nodes
for i=2:N-1
    deltax(i) = x(i) - x(i-1);
    DELTAX(i) = (deltax(i-1)+deltax(i))/2;
end

% Evaluate the initial condition
U_old = 0*ones(N,1);
%% Main Solver
U = U_old; % let the guess of the first solution be the initial condition

% Initialise an array to store the solutons at the times specified by
% plot_times
TEMP_array = zeros(N,length(plot_times));
% Store the initial condition in the array if desired
if plot_times(1) == 0
    TEMP_array(:,1) = U; U_count = 2;
else
    U_count = 1;
end

% Begin the solver
if realtime_plot == true
    figure % Start a new figure
end
% Loop until end time is reached
for t = dt:dt:endtime
    % Call fsolve to solve F(phi) = 0
    [U_fsolve,fval,exitflag,output] = fsolve(@(U) FVM_Advection_Diffusion(U,x,t,U_old,dt,BC,parms,deltax,DELTAX),U_old,optimoptions('fsolve','Display','off'));
    FPARMS = {t,U_old,dt,BC,parms,deltax,DELTAX};
%     y = U;
    x = x;
    m = 1;
    U = Newton_Solve(U,x,m,FPARMS);
    fval =1;
    exitflag = 1;
    output = 1;
    fprintf('t=%g fval=%g exitflag=%g\n',t,norm(fval),exitflag);
%     if exitflag 0, -1, -2, -3 failure to converge solver
    if (exitflag ~= 1)
        error(output)
    end
    U_actual = FVM_Analytic(x,t,parms);
    if realtime_plot == true % Is a plot being produced in realtime?
        hold on
        scatter(x,U,[],[0 0 1]) % plot the numerical solution at each node
        plot(x,U_fsolve)
        hold off
        ylim([0 U0]) % Set the yaxis bounds to 0,1
        xlabel('x') % x-axis label
        ylabel('U(x,t)') % y-axis label
       
        title(['Numerical Solution Time = ' num2str(t)]) % Update the title with the current time
        drawnow % draw the plot in real time
        pause(1e-16) % Plot is updated every 1e-16 seconds or whenever a new solution is ready
% release current plot so the figure is refreshed next iteration
    end
    % Check if the current time is one that is being saved
    if ismembertol(t,plot_times,1e-5)  % check if t == an element in plot_times
        TEMP_array(:,U_count) = U; % save the solution in the array PHI
        U_count = U_count + 1; % increase count for next saved solution
    end
    
    U_old = U; % Update the solution at the previous timestep
end

% Begin plot
figure
hold on
for j=1:length(plot_times) % Loop through plot times
    scatter(x,TEMP_array(:,j),[],rand(1,3)) % Plot Numerical solution
end

% Plot title and axis labels
title('Numerical Solution')
xlim([0 L])
ylabel('U(x,t)')
xlabel('x')

%% Function for evaluating the discretisation
function F = FVM_Advection_Diffusion(U,x,t,U_old,dt,Bcond,parms,deltax,DELTAX)
% Function F that evaluates the vertex-centered finite volume
% discretisation for the linear PDE:
% dU/dt VdU/dx - D*d^2U/dx^2 +lambda*U = 0, 0<x<L, t > 0,
% subject to
% phi(x,0) = f(x) , 0<x<L, initially and
% A*T - B*dT/dx = C(t) on x=0, t>0
% A*T + B*dT/dx = C(t) on x=L, t>0 on the left and right boundaries
% INPUTS
% T - vector of the current estimate for T evaluated at each mesh point
% x - vector of x coordinates for each mesh point
% t - current time
% T_old - vector of the function T evaluated at the previous timestep
% dt - current timestep size
% Bcond - a 2x3 array that stores A,B,C for both boundaries in the form
% above
% parms - A vector containing the temporal weighting scheme, Domain
% length, Number of nodes on the mesh, the diffusivity constant and any
% other parameters [theta L N D V lamda U0 sigma];


% Retrieve the temporal weighting, Domain length, node number and
% diffusivity from parms
theta = parms(1); L = parms(2); N = parms(3); D = parms(4); V = parms(5);
lambda = parms(6); U0 = parms(7); sigma = parms(8);

J = zeros(1,N+1); % vector to store fluxs on each face
Jold = zeros(1,N+1); % vector to store fluxs at previous timestep
F = zeros(N,1); % Vector of values to evaluate function

% Internal Fluxs
for i = 2:N
    J(i) = V*(U(i-1)+ (sigma/2)*(U(i) - U(i-1))) - (D*(U(i)-U(i-1))/deltax(i-1));
    Jold(i) =  V*(U_old(i-1) + (sigma/2)*(U_old(i) - U_old(i-1))) - (D*(U_old(i)-U_old(i-1))/deltax(i-1));
end

% Boundary Fluxs
% Left Boundary x=0
% Have assumed that it is Dirichlet (you can do more complicated conditions
% if you want. 
A = Bcond(1,1); B = Bcond(1,2); C = Bcond(1,3);
if B == 0 % Dirichlet Condition
    F(1) = A*U(1) - C;
else % Neumann or Robin boundary condition
error('Dirichlet Boundary Condition has been assumed')
end

% Right Boundary x = L
% We have assumed that the boundary conditions is dcdx = 0;
    J(N+1) = V*U(N); % Do we need to upwind this term?????
    Jold(N+1) = V*U_old(N);
    % Evaluate the function on the right boundary
    F(N) = U(N) - U_old(N) + dt/DELTAX(N)*(theta*( J(N+1) - J(N)) ...
                           + (1-theta)*( Jold(N+1) - Jold(N))) + dt*lambda*(theta*U(i) + (1-theta)*U_old(i));

% Evaluate the function on the internal nodes
for i = 2:N-1
    F(i) = U(i) - U_old(i) + dt/DELTAX(i)*(theta*( J(i+1) - J(i) + lambda * DELTAX(i)*U(i)) ...
                           + (1-theta)*( Jold(i+1) - Jold(i))  + lambda * DELTAX(i)*U_old(i));
end

end

function U_actual = FVM_Analytic(x,t,parms)
theta = parms(1);
L = parms(2);
N = parms(3);
D = parms(4);
V = parms(5);
lambda = parms(6);
U0 = parms(7);
sigma = parms(8);

U_actual = (U0/2)*(exp((V*(1-lambda).*x)./(2*D)).*...
    erfc((x-V*lambda*t)./(2*sqrt(D*t)) ) +...
    exp((V*(1+lambda).*x)/(2*D)) .* ...
    erfc((x+V*lambda*t)./(2*sqrt(D*t))) ) ;


end
 
%% Function for newton solve
function y = Newton_Solve(y,x,m,FPARMS)

MaxIters = 10; % This is not the max iterations before updating the Jacobian.
func = @(y) Ffunc(y,x,FPARMS);
F = func(y);
k = 1; err = Inf; errold=norm(F,2);
errvecS    = zeros(MaxIters,1);
errvecS(1) = errold;
% tol=tola + tolr*errvecS(1); %Scale Relative tolerance by initial error.

tol = 10^(-8);
rho=1;
h = FPARMS{6};
J = Jacobian(func,y,F);
% Newton Iteration
% fprintf('%2.0f: ||F(y)|| = %1.8e, norm(y-ye,2) = %1.2e \n',k,norm(F,2),norm(y-exact,2));
% Iterate our initial guess y until it satisfies F(y) = 0
while err > tol && k < MaxIters
    if mod(k,m)== 0 || rho>0.5
    J = Jacobian(func,y,F); % This function generates the jacobian. You will have to code this yourself.
    disp("updated the jacobian!")
    end
    dy = J\(-F);
    y = y + dy;
    F = Ffunc(y,x,FPARMS);
    err = norm(F,2);
    k = k + 1;
    errvecS(k) = err; % Store error at each Newton Iteration for plotting purposes. 
    rho = errvecS(k)/errvecS(k-1);
%     fprintf('%2.0f: ||F(y)|| = %1.8e, norm(y-ye,2) = %1.2e \n',k,norm(F,2),norm(y-exact,2));
  
    if err>tol && k == MaxIters
        error('Method has failed as it exceeded the maximum number of iterations and tolerance was not met. \n MaxIters = %2.0f',k)
    end
end
end

function F = Ffunc(y,x,FPARMS)
U = y;
x = x;
t = FPARMS{1};
U_old = FPARMS{2};
dt = FPARMS{3};
Bcond= FPARMS{4};
parms = FPARMS{5};
deltax = FPARMS{6};
DELTAX = FPARMS{7};

F = FVM_Advection_Diffusion(U,x,t,U_old,dt,Bcond,parms,deltax,DELTAX);
end

function J = Jacobian(Ffunc,x,Fx0)
%JACOBIAN Returns the analytical Jacobian of F(x) evaluated at x
% 
% N = size(F,1);
% J = sparse(N,N);
% E = eye(N);
% J(1,1) = 1;
% 
% 
% for j = 2:N-1
%     for i = 2:length(x)-1
%     J(i,j) = (F(i+1)*(E(j)) - F(i))/0.02;
%     end
% end
% J(N,N) = 0;

%%%%%%%%
N = size(x,1);
J = zeros(N,N);
if norm(x,2) == 0
    h = sqrt(eps);
else
    h = sqrt(eps)*norm(x,2);
end

e = zeros(N,1);
for j = 1:N
    e(j) = 1.0;
    J(:,j) = (Ffunc(x+h*e) - Fx0)/h;
    e(j) = 0.0;
end
   
end