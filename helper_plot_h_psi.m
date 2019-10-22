function helper_plot_h_psi(Nz,Nx,figm,X,Z,psi_now,h)

 MAT_h = helper_row2mat(Nz,Nx,h) ;
    MAT_psi = helper_row2mat(Nz,Nx,psi_now) ;
    helper_plotcmap(X,Z,MAT_psi,MAT_h,figm);

end