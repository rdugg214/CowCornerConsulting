function k_out = k_face(k,index,offset)

k_out = ((k(index+offset) + k(index))/2);
end