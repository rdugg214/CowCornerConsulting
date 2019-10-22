function [x,m] = IanGMRES(A,b,x0,M,tol,MaxKrylov,diagnostics)

N = size(A,1);
H = zeros(MaxKrylov+1,MaxKrylov);
V = zeros(N,MaxKrylov+1);

r = b-A*x0;
beta = norm(r,2);
V(:,1) = r/beta;
m = 0;
rnorm = inf;
while rnorm>beta * tol && m<=MaxKrylov
    m = m+1;
    
    V(:,m+1) = A*(M\V(:,m));
    for j = 1:m
        H(j,m) = V(:,j)'*V(:,m+1);
        V(:,m+1) = V(:,m+1) - H(j,m) * V(:,j);
    end
    H(m+1,m) = norm(V(:,m+1),2);
    
    if abs(H(m+1,m))< 1e-14
       fprintf('Invariant Krylov Subspace detected at m=%g\n',m);
       y = H(1:m,1:m)\([beta;zeros(m-1,1)]); %invariant space detected
       break;
    else
        V(:,m+1) = V(:,m+1)/H(m+1,m);
    end
   
    rhs = [beta;zeros(m,1)];
    y = H(1:m+1,1:m)\rhs;
    rnorm = norm(rhs - H(1:m+1,1:m)*y);
    if diagnostics, fprintf('m=%g ||r_m||=%g tol=%g\n',m,rnorm,beta*tol);end
    
end

x = x0+M\(V(:,1:m)*y);
end
