function helper_plotcmap(X,Z,MAT_psi,MAT_h,m)
figure(m)
    subplot(1,4,1)
    surf(X,Z,MAT_h)
%     caxis([-10 2])
    axis([0 500 0 100])
colorbar
    colormap((parula))
    shading interp
    view(2)
   

   subplot(1,4,2)
    surf(X,Z,MAT_psi)
    colormap((parula))
    shading interp
    view(2)
    axis([0 500 0 100])
   caxis([0 1])
  
 
colorbar






end