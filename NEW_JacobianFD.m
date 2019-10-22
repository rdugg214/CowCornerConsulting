function J = NEW_JacobianFD(F, x, Fx0)
% Why recalculate Fx0 all the time should be known from newton iterate.
% ndiag_Jacobian - created by Eamon Conway
% Efficiently calculates the jacobian for a banded system. 
% Future Updates: Probably going to account for variable inputs. 
% If bandwidth not specified calculate it and output as a variable, else,
% allow the code to run with specified band width
% Should there be a flag? I dunno, maybe, not sure when i should throw an 
% error though. 
N = size(x,1);
J = sparse(zeros(N));
ep = sqrt(eps); 
nrmx = norm(x,2);

% Define the size of the shift.
if nrmx == 0
    h = ep;
else
    h = ep*nrmx;
end
I = sparse(eye(N,N));

for j = 1:N
     Fpert = F(x +  (h / norm(I(:,j),2) )*I(:,j));
    Jcol = (Fpert - Fx0)/(h / norm(I(:,j),2) );
    %Update the vector of Diagonals. % Yes, you can vectorise this in a
    %smart way, my tests indicate that the vectorised version is slower. 
    J(:,j) = Jcol;
    
end
end