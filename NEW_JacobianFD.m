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

% Define the size of the shift.
if nrmx == 0
    h = ep;
else
    h = ep*nrmx;
end

FVM_PC = @(h) FVM_pre_calcs(h, dt, t, params);
[k_old, psi_old, Q_old, ~] = FVM_PC(h_old);
F_in = @(h, index, k_new, psi_new, Q, Kzz) FVM_index(h,dt, t, t_old,h_old, params, index, k_old, psi_old, Q_old, k_new, psi_new, Q, Kzz);
for j = 1:N
    x_current = x;
    x_current(j) = x_current(j) + h;
    Fpert = F(x_current);
    Jcol = (Fpert - Fx0)/h;
    J(:,j) = Jcol;
    

%     [k_new, psi_new, Q, Kzz] = FVM_PC(x_current);
%     J(j, j) = (F_in(x_current, [j-params{2}, j-1, j], k_new, psi_new, Q, Kzz) - Fx0(j))/h;
%     diag_offsets = [1, params{2}];
%     for k = 1:length(diag_offsets)
%         neg_index = j - diag_offsets(k);
%         vert_neg_index = [neg_index - diag_offsets(2), neg_index - diag_offsets(1), neg_index];
%         pos_index = j + diag_offsets(k);
%         vert_pos_index = [pos_index - diag_offsets(2), pos_index - diag_offsets(1), pos_index];
%         if (neg_index) > 0
%             J(neg_index, j) = (F_in(x_current, vert_neg_index, k_new, psi_new, Q, Kzz) - Fx0(neg_index))/h;
%         end
%         if (vert_pos_index(3)) <= N
%             J(pos_index, j) = (F_in(x_current, vert_pos_index, k_new, psi_new, Q, Kzz) - Fx0(pos_index))/h;
%         end
%     end
%     Jcol = sparse(Jcol);
%     if ~isequal(Jcol, J(:,j))
%         disp(Jcol);
%         disp(J_test(:,j));
%         error("Jacobian created incorrectly")
%     end
end
end