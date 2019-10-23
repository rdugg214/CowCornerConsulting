function helper_plot_h_psi_av(Nz,Nx,figm,X,Z,psi_now,h,t_hist,psi_av_hist,psi_guess_hist,Ballarr)


 MAT_h = helper_row2mat(Nz,Nx,h) ;
    MAT_psi = helper_row2mat(Nz,Nx,psi_now) ;

    helper_plotcmap(X,Z,MAT_psi,MAT_h,figm);


    subplot(1,4,3)
   
    hold on
    plot(t_hist,psi_av_hist,'Color','r');
     plot(t_hist,psi_guess_hist,'Color','b');
     
     
      subplot(1,4,4)
      plot(t_hist,Ballarr(3,:));
%      'LineWidth',2,
%       ylim([0.04 0.06])
    hold off
    


end