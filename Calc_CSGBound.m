function qe_CSG = Calc_CSGBound(z,hp,Kxx)

DELTACSG = 5000;
HCSG = 20;
qe_CSG = Kxx*(HCSG - (hp+z))/DELTACSG;
end