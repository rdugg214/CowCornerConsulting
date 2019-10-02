% Reset the workspace
clear all, close all, clc;

% This week we are using Newton's method to solve a boundary value problem
fprintf('\n QUESTIONS \n Why is Fsolve bad? Is it slow? \n How has doing our own code improved the total runtime? \n What do you notice about norm(y-ye) for each of the different methods? \n')
% Generate mesh and exact solution
N = 51; a = 1; b=3; h=(b-a)/(N-1);
x=[a:h:b]';
exact=x.^2 + 16./ x;

% Set the tolerances for our newton method
tola=1e-10; tolr=1e-12; %tola is absolute tolerance (used when the solution is close to zero)  and tolr is our relative tolerance.

%% Open Figure 1, This will be used as our error check in all future plots
figure(1)

%% Newton's Method (Using Fsolve)
% In this section is it important to notice that there is no for loop for
% the timestep. This is because after applying the finite volme method to
% the boundary value problem, we arrive at a system of algebraic equations
% and not a system of ordinary differential equations. Therefore, we merely
% have to use fsolve once. Furthermore, the interpretation of y0 varies
% from being the previous time-steps solution to being an initial guess for
% our steady state solution. Note that the initial guess that we have used
% here is the zero vector. What problems could this cause using Newton's
% method? And how could we be smarter? Also, is it worth being smarter? 

fprintf('\n Newton Method with FSolve:\n');
y0 = zeros(N,1); % Define the initial guess y = 0 over the whole domain.
[y,fval,exitflag,output] = fsolve(@(y) Ffunc(y,x,h),y0,optimoptions('fsolve','Display','iter'));

    if (exitflag ~= 1)
        error(output) % Check if fsolve converged to an appropriate solution.
    end

    F = Ffunc(y,x,h);
fprintf('Fsolve: ||F(y)|| = %1.8e, norm(y-ye,2) = %1.2e \n',norm(F,2),norm(y-exact,2));
%% Plot exact vs numerical solution.
figure
plot(x,y,'ko')
hold on
plot(x,exact,'-')
title('Numerical vs Exact Solution (Fsolve)')
hold off

%% Full Newton's Method 
MaxIters = 10;
fprintf('\n (Full) Newton Method:\n\n');
y = zeros(N,1); % Can we have a smarter initial guess?
F = Ffunc(y, x, h);
errvecN = zeros(MaxIters,1);
errvecN(1) = norm(F,2);
tol=tola + tolr*errvecN(1); %Scale Relative tolerance by initial error. 
k = 1; err = Inf;
fprintf('%2.0f: ||F(y)|| = %1.8e, norm(y-ye,2) = %1.2e \n',k,norm(F,2),norm(y-exact,2));
% Iterate our initial guess y until it satisfies F(y) = 0
while err > tol && k < MaxIters
    J = Jacobian(y,h); % This function generates the jacobian. You will have to code this yourself. 
    dy = J\(-F);
    y = y + dy;
    F = Ffunc(y,x,h);
    err = norm(F,2);
    k = k + 1;
    errvecN(k) = err; % Store error at each Newton Iteration for plotting purposes. 
    fprintf('%2.0f: ||F(y)|| = %1.8e, norm(y-ye,2) = %1.2e \n',k,norm(F,2),norm(y-exact,2));
    
    if err>tol && k == MaxIters
        error('Newton Method has failed as it exceeded the maximum number of iterations and tolerance was not met. \n MaxIters = %2.0f',k)
    end
end

%% Plot exact vs numerical solution.
figure
plot(x,y,'ko')
hold on
plot(x,exact,'-')
title('Numerical vs Exact Solution (Full)')
hold off

%% Plot convergence of Full Newton's Method
figure(1) % Remember that we initialised figure 1 to be the convergence plots. 
logerr = log10(errvecN);
plot(logerr(1:k-1),logerr(2:k),'bx','MarkerSize',8)
hold on
% Line of best fit: log(err^(k+1)) = c(1) + c(2)*log(err^(k))
% Slope c(2) is the convergence rate
A = [ones(k-1,1), logerr(1:k-1)];
b = logerr(2:k);
c =A\b;
NewtonRate = c(2);
ll=@(x) c(1)+c(2)*x;
plot(logerr(1:k-1),ll(logerr(1:k-1)),'b-');
hold off

%% Chord Method
fprintf('\n Chord Method:\n\n');
MaxIters = 200; % Max number of newton steps.
y = zeros(N,1);
F = Ffunc(y, x, h);
J = Jacobian(y, h); % This function generates the jacobian. You will have to code this yourself. 
[L,U] = lu(J); %Why would you want to do this for the Chord method, but not the Full Newon Method? Could this be extended to timestepping....? 
errvecC = zeros(MaxIters,1);
errvecC(1) = norm(F,2);
k = 1; err = Inf; 
tol=tola + tolr*errvecC(1); %Scale Relative tolerance by initial error.
fprintf('%2.0f: ||F(y)|| = %1.8e, norm(y-ye,2) = %1.2e \n',k,norm(F,2),norm(y-exact,2));
% Newton Iteration
while err > tol && k < MaxIters
    dy = L\(-F); % Forward sub
    dy = U\dy; % Backward sub
    y = y + dy;
    F = Ffunc(y,x,h);
    err = norm(F,2); 
    k = k + 1;
    errvecC(k) =  err; % Store error at each Newton Iteration for plotting purposes. 
    fprintf('%2.0f: ||F(y)|| = %1.8e, norm(y-ye,2) = %1.2e \n',k,norm(F,2),norm(y-exact,2));
    if err>tol && k == MaxIters
        error('Chord Method has failed as it exceeded the maximum number of iterations and tolerance was not met. \n MaxIters = %2.0f',k)
    end
    
end


%% Plot exact vs numerical solution.
figure
plot(x,y,'ko')
hold on
plot(x,exact,'-')
title('Numerical vs Exact Solution (Chord)')
hold off

%% Plot convergence of Chord Method
figure(1)
hold on
logerr = log10(errvecC);
plot(logerr(1:k-1),logerr(2:k),'ro')
% Line of best fit: log(err^(k+1)) = c(1) + c(2)*log(err^(k))
% Slope c(2) is the convergence rate
A = [ones(k-1,1), logerr(1:k-1)];
b = logerr(2:k);
c = A\b;
ChordRate = c(2);
ll=@(x) c(1)+c(2)*x;
plot(logerr(1:k-1),ll(logerr(1:k-1)),'r-');
hold off

%% Newton-Shamanskii Method
% Set relative and Absolute Tolerances
MaxIters = 10; % This is not the max iterations before updating the Jacobian.
m = 3; % Iterations before updating Jacobian.
fprintf('\n Shamanskii Method (m = %i):\n',m);
k = 1; err = Inf; errold=norm(F,2);
errvecS    = zeros(MaxIters,1);
errvecS(1) = errold;
tol=tola + tolr*errvecS(1); %Scale Relative tolerance by initial error.
rho=1;

% Newton Iteration



%% Plot exact vs numerical solution.

%% Plot convergence of Newton-Shamanskii Method
figure(1)
hold on

% Line of best fit: log(err^(k+1)) = c(1) + c(2)*log(err^(k))
% Slope c(2) is the convergence rate

ShamanskiiRate = nan;

% Plot straight line showing convergence rate. 

hold off

%% Inexact Newton
fprintf('\n Inexact Newton Method:\n');
y = zeros(N,1);
F = Ffunc(y, x, h);
MaxIters = 10; % Maximum newton steps.
errvecIN    = zeros(MaxIters,1);
errvecIN(1) = norm(F,2);
k = 1; err = Inf;
tol=tola + tolr*errvecIN(1);%Scale Relative tolerance by initial error.

% Newton Iteration

%% Plot exact vs numerical solution.

%% Plot convergence of Inexact Newton Method
figure(1)
hold on

% Line of best fit: log(err^(k+1)) = c(1) + c(2)*log(err^(k))
% Slope c(2) is the convergence rate

InexactNewtonRate = nan;

% Plot straight line showing convergence rate. 

hold off

%% Label the axis of our convergence plots.
xlabel('$\textrm{log}_{10}(\|F(\mathbf{x}^{(k)})\|_{2})$',...
    'Interpreter','LaTeX')
ylabel('$\textrm{log}_{10}(\|F(\mathbf{x}^{(k+1)})\|_{2})$',...
    'Interpreter','LaTeX')
set(gca,'FontSize',14)
h = legend('(Full) Newton','Chord','Shamanskii','Inexact Newton',...
    'Location','NorthWest');
set(h,'FontSize',14)


%% Convergence Rates
fprintf('\nConvergence Rates:\n')
fprintf('(Full) Newton Method:         %1.3f\n',NewtonRate);
fprintf('Chord Method:                 %1.3f\n',ChordRate);
fprintf('Shamanskii Method (m = %i):    %1.3f\n',m,ShamanskiiRate);
fprintf('Inexact (Full) Newton Method: %1.3f\n',InexactNewtonRate);

