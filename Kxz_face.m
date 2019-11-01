function K_OUT = Kxz_face(Kxz,offset,index,boundary)
%     
% if index == 229
%     disp('test')
% end

% if ~boundary 
% % K_OUT=Kxz(index+offset);
% K_OUT = Kxz(index);

% else
K_OUT = ((Kxz(index+offset)+Kxz(index)))/2;
% end
end