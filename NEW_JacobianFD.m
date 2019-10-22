function J = NEW_JacobianFD(F, x, Fx0, F_in)
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

J_test = sparse(zeros(N));

% Define the size of the shift.
if nrmx == 0
    h = ep;
else
    h = ep*nrmx;
end
I = sparse(eye(N,N));

for j = 1:N
%      Fpert = F(x +  (h / norm(I(:,j),2) )*I(:,j));
     Fpert = F(x);
%      Fpert = F(x +  h*I(:,j));
    Jcol = (Fpert - Fx0)/(h / norm(I(:,j),2) );
    %Update the vector of Diagonals. % Yes, you can vectorise this in a
    %smart way, my tests indicate that the vectorised version is slower. 
    J(:,j) = Jcol;
    
    
    J_test(j, j) = F_in(x, j);
%     J_test(j, j) = F_in(x + h*I(:,j), j);
%     J_test(j, j) = (F_in(x + h*I(:,j), j) - Fx0(j))/h;
    diag_offsets = [1, 11];
    for k = 1:length(diag_offsets)
        if (j - diag_offsets(k)) > 0
            J_test(j - diag_offsets(k), j) = F_in(x, j - diag_offsets(k));
%             J_test(j - diag_offsets(k), j) = (F_in(x, j - diag_offsets(k)) - Fx0(j))/h;
        end
        if (j + diag_offsets(k)) <= N
            J_test(j + diag_offsets(k), j) = F_in(x, j + diag_offsets(k));
%             J_test(j + diag_offsets(k), j) = (F_in(x, j + diag_offsets(k)) - Fx0(j))/h;
        end
    end
    if ~isequal(Fpert, J_test(:,j))
        disp(J(:,j));
        disp(J_test(:,j));
        error("This is where is should fucking stop")
    end
end
disp(J)
figure
spy(J)
figure
spy(J - J_test)

end