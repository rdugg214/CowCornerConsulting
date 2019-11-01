function [Q,nKzz] = Calc_Q(h,x,z,dt,psi,psi_sat,t,t_on,Kzz,DX,DZ,simple,Pr,prediction_data)

Q = zeros(size(h));


% % Evapotrans
L2 = 100;
if simple
Q_evap = @(Ri,li,Z)  -(((Ri/100)*(1 +  cos((2*pi*t)/365 ))*1.803 * 0.001 )  *(Z-L2+li).^2)/(li^2);
else
Q_evap = @(Ri,li,Z)  -(((Ri/100)*rainfall(1,t,prediction_data)) *(Z-L2+li).^2)/(li^2);
end
 loc =(50<=x & x <= 100) & (85 <= z & z<= 100) & (psi./psi_sat > 0.5) & ~simple;
Q(loc) = Q_evap(5,15,z(loc));

loc = (100<=x & x <= 300) & (95 <= z & z<= 100) & (psi./psi_sat > 0.5) & ~simple;
Q(loc) = Q_evap(2.5,5,z(loc));

loc = (300<=x & x <= 500) & (90 <= z & z<= 100) & (psi./psi_sat > 0.5) & ~simple;
Q(loc) = Q_evap(3.5,10,z(loc));


% Pumping
nKzz = Kzz;
if simple 
Tarf =0.6544;
else
Tarf = prediction_data(ceil(t/365))/1000;
end

loc = x==100  & 55 <= z & z<= 75 & t >= t_on ;%& (psi./psi_sat)>0.11;
fr = 50;
Q(loc) = -(Tarf * Pr)./((365)*sum(loc));

nKzz(loc) = Kzz(loc) *fr;

end