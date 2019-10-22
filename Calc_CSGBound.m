function qe_CSG = Calc_CSGBound(z,hp,Kxx,t,t_on_CSG)
qe_CSG = 0;
if z <=5 && t >= t_on_CSG
DELTACSG = 5000;
HCSG = 20;
qe_CSG = (Kxx*(HCSG - (hp+z))/DELTACSG);
end
end