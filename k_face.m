function k_out = k_face(k,d,index,offset)

sigma = 0;
k_out = ((k(index+offset) + k(index))/2); %arithmetic average
% i = index+offset;
% if offset == 1
%     k_out = k(index) - (sigma/2)*(k(index) - k(index+1));
% elseif offset == -1
%      k_out = k(index-1) - (sigma/2)*(k(index-1) - k(index));
% elseif offset <0 %-NX
%      k_out = k(index) - (sigma/2)*(k(index) - k(index+offset));
% else  %NX
%      k_out = k(index+offset) - (sigma/2)*(k(index+offset) - k(index));
% end
% if offset>0
% k_out =     d(index) * ((d(index)/2)/k(index) + (d(index)/2)/k(i) )^(-1);
% else
% k_out = d(i) * ((d(i)/2)/k(index) + (d(i)/2)/k(i) )^(-1); %harmonic
% end

end