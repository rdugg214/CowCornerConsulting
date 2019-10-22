function output = GP_sym(start,endPoint,n,r)
% Making sure n is not a decimal or an even number. Adjusting if so. 
    l = endPoint-start; % length
    n= round(n);
    if (mod(n,2) == 0)
        if (mod(n,2) == 0)
            n = n+1;
            fprintf('Adjusted number of points to n = %d\n', n)
        end
    end
% Calculate delta x for each node
    dx = (l/2)*(1-r) / (1- r^((n+1)/2-1));
    output = [start];
% Generate symmetric geometric progression. 
    for i  = 2:n
        if i <= (n+1)/2 
         output(i) = output(i-1) + dx * r^(i-2);
        elseif i >= (n+1)/2 +1
         output(i) = output(i-1) + dx * r^(n-i);
        end       
    end
end