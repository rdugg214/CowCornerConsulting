function [x_out,success] = NEW_Newton_Solver(F,x0,Jacobian, type)
%% Inexact Newton
success = true;
tola = 1e-6;
tolr = 1e-6;
fprintf('\n Inexact Newton Method:\n');
x = x0;
Resid = F(x);
MaxIters = 10;
tol = tola + norm(x)*tolr;
% Newton Iteration
m = 3; % Number of steps until you update jacobian. 
preverr = 0;
maxerr = 1e0;
err =Inf;
if (type ~= "Full Newton")
    J = Jacobian(F,x,Resid);  
end
Method = "GMRES Test";
k = 0;
while err > tol && k < MaxIters
    if CalculateJacobian(type, k, m, err, preverr)
        J = Jacobian(F,x,Resid); 
        M = Find_Precon_Mat(J, "ILU", 5);
    end
    if Method == "Normal"
        dx = J\(-Resid);
    else 
    
        dx = IanGMRES(J,-Resid,x,M,tol,20,0);
    end
    lambda = 1;
    alpha = 1e-2;
    xt = x + lambda*dx;
%     line_search_count = 0;
%     F_norm = norm(F(x),1)^2;
%     Fxt_norm = norm(F(xt),1)^2;
%     while Fxt_norm >  (1-2*alpha*lambda)*F_norm
%         lambda_star = (F_norm * lambda^2)/(Fxt_norm + F_norm * (2 * lambda - 1));
%         if lambda_star < lambda * 0.1
%             lambda = lambda * 0.1;
%         elseif lambda_star > lambda * 0.5
%             lambda = lambda * 0.5;
%         else
%             lambda = lambda_star;
%         end
%         xt = x+lambda*dx;
%         Fxt_norm = norm(F(xt),1)^2;
%         fprintf('Norm xy %1.8e --- Norm x %1.8e \n', norm(F(xt),1)^2,  (1-2*alpha*lambda)*norm(F(x),1)^2);
%         line_search_count = line_search_count + 1;
%     end
%     fprintf('Line search took %d iterations to complete \n', line_search_count);
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