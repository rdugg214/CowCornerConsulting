function J = NEW_JacobianFD(F, x, Fx0, params, dt, t, t_old, h_old)
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

old_method_time = 0;
new_method_time = 0;

% FVM_PC = @(h) FVM_pre_calcs(h, dt, t, params);
% [k_old, psi_old, Q_old, ~] = FVM_PC(h_old);
% F_in = @(h, index, k_new, psi_new, Q, Kzz) FVM_index(h,dt, t, t_old,h_old, params, index, k_old, psi_old, Q_old, k_new, psi_new, Q, Kzz);
for j = 1:N
%     tic;
    Fpert = F(x +  h*I(:,j));
    Jcol = (Fpert - Fx0)/(h / norm(I(:,j),2) );
    J(:,j) = Jcol;
%     old_method_time = old_method_time + toc;
    
%     tic;
%     [k_new, psi_new, Q, Kzz] = FVM_PC(x + h*I(:,j));
%     J_test(j, j) = (F_in(x + h*I(:,j), j, k_new, psi_new, Q, Kzz) - Fx0(j))/h;
%     diag_offsets = [1, params{2}];
%     for k = 1:length(diag_offsets)
%         vert_neg_index = j - diag_offsets(k);
%         vert_pos_index = j + diag_offsets(k);
%         if (vert_neg_index) > 0
%             J_test(vert_neg_index, j) = (F_in(x + h*I(:,j), vert_neg_index, k_new, psi_new, Q, Kzz) - Fx0(vert_neg_index))/h;
%         end
%         if (vert_pos_index) <= N
%             J_test(vert_pos_index, j) = (F_in(x + h*I(:,j), vert_pos_index, k_new, psi_new, Q, Kzz) - Fx0(vert_pos_index))/h;
%         end
%     end
   
%     new_method_time = new_method_time + toc;
%     Jcol = sparse(Jcol);
%     if ~isequal(Jcol, J_test(:,j))
%         disp(Jcol);
%         disp(J_test(:,j));
%         error("This is where is should fucking stop")
%     end
end
% disp(old_method_time);
% disp(new_method_time);
% spy(J - J_test)
%  J = J_test;
end