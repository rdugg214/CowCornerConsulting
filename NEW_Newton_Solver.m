function [x_out,success] = NEW_Newton_Solver(F,x,Jacobian)
%% Inexact Newton
success = true;
tola = 1e-6;
tolr = 1e-6;
fprintf('\n Inexact Newton Method:\n');
Resid = F(x);
MaxIters = 15;
tol = tola + norm(x)*tolr;
% Newton Iteration
m = 3; % Number of steps until you update jacobian. 
preverr = 0;
maxerr = 1e0;
err =Inf;
J = Jacobian(F,x,Resid);  
k = 0;
setup.type = 'nofill';
setup.milu = 'off';
setup.droptol = 1e-6;
   M = ilu(sparse(J),setup);
while err > tol && k < MaxIters
    if mod(k,m)==0 || preverr < err
    J = Jacobian(F,x,Resid);  
   
    end
    dx = J\(-Resid);
%     dx = IanGMRES(J,-Resid,x,M,tol,20,0);
     lambda = 1;
    alpha = 1e-2;
    x = x + lambda*dx;
    xt = x;
%     while norm(F(xt),1)^2 >  (1-2*alpha*lambda)*norm(F(x),1)^2
%         sigma = (0.5-0.1) * rand + 0.1;
%         lambda = sigma*lambda;
%         xt = x+lambda*dx;
%         fprintf('Norm xy %1.8e --- Norm x %1.8e \n', norm(F(xt),1)^2,  (1-2*alpha*lambda)*norm(F(x),1)^2);
%     end
    x = xt;
    k = k + 1;
    Resid = F(x);
    preverr = err;
    err = norm(Resid,2);
    fprintf('%2.0f: ||F(y)|| = %1.8e \n',k,norm(Resid,2));
    if err >maxerr || isnan(err) || k >= MaxIters || isinf(err)
        success = false;
         break
    end
end

x_out = x;