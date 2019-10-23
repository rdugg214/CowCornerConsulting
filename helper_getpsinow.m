function psi_now = helper_getpsinow(h, alpha,n,m,psi_res,psi_sat,x,z,dx,dz,hetgen)
S_now = CalcS(h, alpha, n, m);
k_now = Calck(h, S_now, m,x,z,dx,dz,hetgen);
psi_now = CalcPsi(h, S_now, psi_res, psi_sat,x,z,dx,dz,hetgen);
end