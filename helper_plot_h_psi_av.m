function helper_plot_h_psi_av(Nz,Nx,figm,X,Z,psi_now,h,t_hist,psi_av_hist,psi_guess_hist,mode)

 if mode == 1
 MAT_h = helper_row2mat(Nz,Nx,h) ;
    MAT_psi = helper_row2mat(Nz,Nx,psi_now) ;
 else
     MAT_h = reshape(h,Nx,Nz)';
     MAT_psi = reshape(psi_now,Nx,Nz)';
     end
    helper_plotcmap(X,Z,MAT_psi,MAT_h,figm);


    subplot(1,3,3)
   
    hold on
    plot(t_hist,psi_av_hist,'Color','r');
     plot(t_hist,psi_guess_hist,'Color','b');
%      'LineWidth',2,
%       ylim([0.04 0.06])
    hold off
    


end