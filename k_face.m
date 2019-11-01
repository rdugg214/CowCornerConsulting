function k_out = k_face(k,d,index,offset)

sigma = 0;
sigmax = 1;
if offset == 1
    k_out = k(index) - (sigmax/2)*(k(index) - k(index+1));
elseif offset == -1
     k_out = k(index-1) - (sigmax/2)*(k(index-1) - k(index));
elseif offset <0 %-NX
     k_out = k(index) - (sigma/2)*(k(index) - k(index+offset));
else  %NX
     k_out = k(index+offset) - (sigma/2)*(k(index+offset) - k(index));
end

end