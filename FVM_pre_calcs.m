function [k, psi, Q, Kzz] = FVM_pre_calcs(h, dt, t, params)
    alpha = params{4};      %vector size h
    n = params{5};          %vector size h
    m = params{6};          %vector size h
    psi_res = params{7};    %vector size h
    psi_sat = params{8};    %vector size h
    x = params{9};          %vector size h
    z = params{10};         %vector size h
    Kzz = params{12};       %vector size h
    dx = params{13};        %vector size h
    dz = params{14};        %vector size h
    DELTAX =  params{15};   %vector size h
    DELTAZ = params{16};    %vector size h
    t_on_PUMP = params{18};
    simple = params{19};
    Pr = params{20};
    hetgen = params{21};
    prediction_data = params{22};
    DELCSG = params{23};
    
    S = CalcS(h, alpha, n, m);
    k = Calck(h, S, m, x, z, dx, dz, hetgen);
    psi = CalcPsi(h, S, psi_res, psi_sat,x,z,dx,dz,hetgen);
    [Q,Kzz]=  Calc_Q(h,x,z,dt,psi,psi_sat,t,t_on_PUMP,Kzz,DELTAX,DELTAZ,simple,Pr,prediction_data);
end