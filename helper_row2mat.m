function MAT = helper_row2mat(Nz,Nx,ROW) 
for i = 1:Nz
    for j = 1:Nx
        index = (i-1)*Nx + j;
        MAT(i,j) = ROW(index);
    end
end
end