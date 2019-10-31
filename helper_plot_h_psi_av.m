function helper_plot_h_psi_av(Nz,Nx,figm,X,Z,psi_now,h,t_hist,psi_av_hist,psi_guess_hist,t,dt,simple,DELCSG,Pr)


MAT_h = helper_row2mat(Nz,Nx,h) ;
MAT_psi = helper_row2mat(Nz,Nx,psi_now) ;

figure(figm)
clf

subplot(1,3,1)
surf(X,Z,MAT_h)
caxis([-10 2])
axis([0 500 0 100])
ylabel('Z (m)')
xlabel('X (m)')
h = colorbar;
colormap((parula))
shading interp
view(2)
ylabel(h,'Pressure head (h) (m)')



subplot(1,3,2)
surf(X,Z,MAT_psi)

ylabel('Z (m)')
xlabel('X (m)')

h = colorbar;
axis([0 500 0 100])
caxis([0 1])
colormap((parula))
shading interp
view(2)
ylabel(h,'Percent Saturation (\psi) (%)')
title(sprintf('t = %.2f (years) (current dt = %.2f (days))',t/365,dt));

subplot(1,3,3)
hold on
plot(t_hist/365,psi_av_hist,'Color','r');
if simple
    plot(t_hist/365,psi_guess_hist,'Color','b');
end
xlabel('time (years)')
ylabel('Average Water Content')
title(sprintf('Parameters DELTACSG = %g, Pumping R = %0.2f',DELCSG,Pr));
end