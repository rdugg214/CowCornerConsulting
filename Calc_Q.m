function [Q,nKzz] = Calc_Q(h,x,z,dt,psi,psi_sat,t,t_on,Kzz,DX,DZ)

Q = zeros(size(h));


% % Evapotrans
L2 = 100;
Q_evap = @(Ri,li,Z)  -(((Ri/100)*(1 +  cos((2*pi*t)/365 ))*1.803 * 0.001 )  *(Z-L2+li).^2)/(li^2);
 loc =(50<=x & x <= 100) & (85 <= z & z<= 100) & (psi./psi_sat > 0.5);
Q(loc) = Q_evap(5,15,z(loc));

loc = (100<=x & x <= 300) & (95 <= z & z<= 100) & (psi./psi_sat > 0.5);
Q(loc) = Q_evap(2.5,5,z(loc));

loc = (300<=x & x <= 500) & (90 <= z & z<= 100) & (psi./psi_sat > 0.5);
Q(loc) = Q_evap(3.5,10,z(loc));


% Pumping
nKzz = Kzz;
Tarf =0.6544;
loc = x==100  & 55 <= z & z<= 75 & t >= t_on ;%& (psi./psi_sat)>0.11;
fr = 50;
Pr = 0.25;
Q(loc) = -(Tarf * Pr)./((365)*sum(loc).*DX(loc).*DZ(loc));
% loc = 80< x & x<120  & 40 <= z  & z<= 90;
nKzz(loc) = Kzz(loc) *fr;
% % % HAVE TO FIX FLUX AT PUMP

end