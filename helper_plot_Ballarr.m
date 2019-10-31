function helper_plot_Ballarr(t_hist,Ballarr,figb,DELCSG,Pr,t,dt)

% 1- riverloc = z >=80 & x ==0 ;
% 2- CSGloc = z <=5 & x ==500;
% 3-rainloc = z ==100;
% 4-evaploc = ((50<=x & x <= 100) & (85 <= z & z<= 100)) |  ((100<=x & x <= 300) & (95 <= z & z<= 100)) ...
%     | ((300<=x & x <= 500) & (90 <= z & z<= 100));
% 5- pumploc = x==100  & 55 <= z & z<= 75;
figure(figb)
clf
t_hist = t_hist/365;
subplot(1,2,1)
hold on
plot(t_hist,-Ballarr(3,:))
plot(t_hist,Ballarr(1,:),'b')
legend('Rain','River')
xlabel('Time (years)')
ylabel('Flux (m/d)')
title(sprintf('t = %.2f (years) (current dt = %.2f (days))',t/365,dt));

subplot(1,2,2)
hold on
plot(t_hist,Ballarr(2,:),'g')
plot(t_hist,Ballarr(4,:),'r')
plot(t_hist,Ballarr(5,:),'k')
hold off
xlabel('Time (years)')
ylabel('Flux (m/d)')
% yyaxis left
% plot(t_hist,Ballarr([1 2 4 5],:))
legend('CSG','Evap','Pump')
title(sprintf('Parameters DELTACSG = %g, Pumping R = %0.2f',DELCSG,Pr));
end