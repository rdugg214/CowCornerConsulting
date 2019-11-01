function qe_CSG = Calc_CSGBound(z,hp,Kxx,t,t_on_CSG,DELCSG)
qe_CSG = 0;
if z <=5 && t >= t_on_CSG  
HCSG = 20;
qe_CSG = -(Kxx*(HCSG - (hp+z))/DELCSG);
end
end