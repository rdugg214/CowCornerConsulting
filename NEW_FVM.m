 function [F] = NEW_FVM(h,dt, t, t_old,h_old, params,node_params)
CSG_HEIGHT = 30;
% Parameters IN
N = params{1};          %value
Nx = params{2};         %value
Nz = params{3};         %value
alpha = params{4};      %vector size h
n = params{5};          %vector size h
m = params{6};          %vector size h
psi_res = params{7};    %vector size h
psi_sat = params{8};    %vector size h
X = params{9};          %vector size h
Z = params{10};         %vector size h
Kxx = params{11};       %vector size h
Kzz = params{12};       %vector size h
dxe=params{13};
dxw=params{14};
dzn =params{15};
dzs=params{16};
DZ=params{17};
DX=params{18};
% node_num = params{19};

IN_P = node_params{1};
IN_N= node_params{2};
IN_E= node_params{3};
IN_W= node_params{4};
IN_S= node_params{5};
P = (1:N)';
TOP_P =node_params{6};
BOT_P = node_params{7};
L_P= node_params{8};
R_P=node_params{9};
node_num = node_params{10};
node_num= node_num(:);

TL_P =node_params{11};
TR_P =node_params{12};
BL_P =node_params{13};
BR_P =node_params{14};
X = X(:);
Z = Z(:);

S = CalcS(h, alpha, n, m);
S_old = CalcS(h_old, alpha, n, m);
k = Calck(h, S, m);
k_old = Calck(h_old, S_old, m);
psi = CalcPsi(h, S, psi_res, psi_sat);
psi_old = CalcPsi(h_old, S_old, psi_res, psi_sat);

Q = -Calc_Q(h,X,Z,dt); %To be updated
Q_old = -Calc_Q(h_old,X,Z,dt);
theta = 0.5;

F = zeros(N,1);
flux_n = F;
flux_e = F;
flux_n_old = F;
flux_e_old = F;

%% Region 1
flux_n(IN_P) = -k(IN_P) .* Kzz(IN_P) .* ((h(IN_P+1) - h(IN_P))./dzn(IN_P) + 1);
flux_e(IN_P) = -k(IN_P) .* Kzz(IN_P) .* ((h(IN_P+Nz) - h(IN_P))./dxe(IN_P));
flux_n_old(IN_P) = -k_old(IN_P) .* Kzz(IN_P) .* ((h_old(IN_P+1) - h_old(IN_P))./dzn(IN_P) + 1);
flux_e_old(IN_P) = -k_old(IN_P) .* Kzz(IN_P) .* ((h_old(IN_P+Nz) - h_old(IN_P))./dxe(IN_P));
F(IN_P) = psi(IN_P) - psi_old(IN_P) + ...
    (theta *dt *( (1./DX(IN_P)) .*(flux_e(IN_P) - flux_e(IN_P-Nz))      +...
    (1./DZ(IN_P)) .* (flux_n(IN_P) - flux_n(IN_P-1)))  +Q(IN_P)) + ...
    ((1-theta) *dt *( (1./DX(IN_P)) .*(flux_e_old(IN_P) - flux_e_old(IN_P-Nz))      +...
    (1./DZ(IN_P)) .* (flux_n_old(IN_P) - flux_n_old(IN_P-1) ) )  +Q_old(IN_P)) ;


%% TOP face
flux_e(TOP_P) = -k(TOP_P) .* Kzz(TOP_P) .* ((h(TOP_P+Nz) - h(TOP_P))./dxe(TOP_P));
flux_e_old(TOP_P) = -k_old(TOP_P) .* Kzz(TOP_P) .* ((h_old(TOP_P+Nz) - h_old(TOP_P))./dxe(TOP_P));

F(TOP_P) = psi(TOP_P) - psi_old(TOP_P) + ...
    (theta *dt *( (1./DX(TOP_P)) .*(flux_e(TOP_P) - flux_e(TOP_P-Nz))      +...
    (1./DZ(TOP_P)) .* (Calc_RainBound(t,h(TOP_P))- flux_n(TOP_P-1))) + Q(TOP_P)) + ...
    ((1-theta) *dt *( (1./DX(TOP_P)) .*(flux_e_old(TOP_P) - flux_e_old(TOP_P-Nz))      +...
    (1./DZ(TOP_P)) .* (Calc_RainBound(t_old,h_old(TOP_P)) - flux_n_old(TOP_P-1) ) )+ Q_old(TOP_P) ) ;

%% BOT Face
flux_n(BOT_P) = -k(BOT_P) .* Kzz(BOT_P) .* ((h(BOT_P+1) - h(BOT_P))./dzn(BOT_P) + 1);
flux_e(BOT_P) = -k(BOT_P) .* Kzz(BOT_P) .* ((h(BOT_P+Nz) - h(BOT_P))./dxe(BOT_P));
flux_n_old(BOT_P) = -k_old(BOT_P) .* Kzz(BOT_P) .* ((h_old(BOT_P+1) - h_old(BOT_P))./dzn(BOT_P) + 1);
flux_e_old(BOT_P) = -k_old(BOT_P) .* Kzz(BOT_P) .* ((h_old(BOT_P+Nz) - h_old(BOT_P))./dxe(BOT_P));
F(BOT_P) = psi(BOT_P) - psi_old(BOT_P) + ...
    (theta *dt *( (1./DX(BOT_P)) .*(flux_e(BOT_P) - flux_e(BOT_P-Nz))      +...
    (1./DZ(BOT_P)) .* (flux_n(BOT_P))) +Q(BOT_P)) + ...
    ((1-theta) *dt *( (1./DX(BOT_P)) .*(flux_e_old(BOT_P) - flux_e_old(BOT_P-Nz))      +...
    (1./DZ(BOT_P)) .* (flux_n_old(BOT_P)) )+ Q_old(BOT_P))  ;

%% EAST Face
flux_n(R_P) = -k(R_P) .* Kzz(R_P) .* ((h(R_P+1) - h(R_P))./dzn(R_P) + 1);
flux_n_old(R_P) = -k_old(R_P) .* Kzz(R_P) .* ((h_old(R_P+1) - h_old(R_P))./dzn(R_P) + 1);
F(R_P) = psi(R_P) - psi_old(R_P) + ...
    (theta *dt *( (1./DX(R_P)) .*(Calc_RightBound(Z(R_P),h(R_P),Kxx(R_P),CSG_HEIGHT) - flux_e(R_P-Nz))      +...
    (1./DZ(R_P)) .* (flux_n(R_P) - flux_n(R_P-1)))+ Q(R_P) ) + ...
    ((1-theta) *dt *( (1./DX(R_P)) .*(Calc_RightBound(Z(R_P),h_old(R_P),Kxx(R_P),CSG_HEIGHT) - flux_e_old(R_P-Nz))      +...
    (1./DZ(R_P)) .* (flux_n_old(R_P) - flux_n_old(R_P-1) ) )+Q_old(R_P) ) ;

%% LEFT 
flux_n(L_P) = -k(L_P) .* Kzz(L_P) .* ((h(L_P+1) - h(L_P))./dzn(L_P) + 1);
flux_e(L_P) = -k(L_P) .* Kzz(L_P) .* ((h(L_P+Nz) - h(L_P))./dxe(L_P));
flux_n_old(L_P) = -k_old(L_P) .* Kzz(L_P) .* ((h_old(L_P+1) - h_old(L_P))./dzn(L_P) + 1);
flux_e_old(L_P) = -k_old(L_P) .* Kzz(L_P) .* ((h_old(L_P+Nz) - h_old(L_P))./dxe(L_P));
F(L_P) = psi(L_P) - psi_old(L_P) + ...
    (theta *dt *( (1./DX(L_P)) .*(flux_e(L_P) - Calc_LeftBound(Z(L_P),h(L_P)))      +...
    (1./DZ(L_P)) .* (flux_n(L_P) - flux_n(L_P-1))) +Q(L_P)) + ...
    ((1-theta) *dt *( (1./DX(L_P)) .*(flux_e_old(L_P) - Calc_LeftBound(Z(L_P),h_old(L_P)))      +...
    (1./DZ(L_P)) .* (flux_n_old(L_P) - flux_n_old(L_P-1) ) ) +Q_old(L_P)) ;

%% TOP LEFT
flux_e(TL_P) = -k(TL_P) .* Kzz(TL_P) .* ((h(TL_P+Nz) - h(TL_P))./dxe(TL_P));
flux_e_old(TL_P) = -k_old(TL_P) .* Kzz(TL_P) .* ((h_old(TL_P+Nz) - h_old(TL_P))./dxe(TL_P));

F(TL_P) = psi(TL_P) - psi_old(TL_P) + ...
    (theta *dt *( (1./DX(TL_P)) .*(flux_e(TL_P) -  Calc_LeftBound(Z(TL_P),h(TL_P))      +...
    (1./DZ(TL_P)) .* (Calc_RainBound(t,h(TL_P)) - flux_n(TL_P-1))))  +Q(TL_P)) + ...
    ((1-theta) *dt *( (1./DX(TL_P)) .*(flux_e_old(TL_P) -  Calc_LeftBound(Z(TL_P),h_old(TL_P))      +...
    (1./DZ(TL_P)) .* (Calc_RainBound(t_old,h_old(TL_P)) - flux_n_old(TL_P-1) ) ))  +Q_old(TL_P)) ;
% 
%% TOP RIGHT
F(TR_P) = psi(TR_P) - psi_old(TR_P) + ...
    (theta *dt *( (1./DX(TR_P)) .*(- flux_e(TR_P-Nz))      +...
    (1./DZ(TR_P)) .* (Calc_RainBound(t,h(TR_P))- flux_n(TR_P-1)))  +Q(TR_P)) + ...
    ((1-theta) *dt *( (1./DX(TR_P)) .*( - flux_e_old(TR_P-Nz))      +...
    (1./DZ(TR_P)) .* (Calc_RainBound(t_old,h_old(TR_P)) - flux_n_old(TR_P-1) ) ) +Q_old(TR_P)) ;
% 
%% BOT LEFT
flux_n(BL_P) = -k(BL_P) .* Kzz(BL_P) .* ((h(BL_P+1) - h(BL_P))./dzn(BL_P) + 1);
flux_e(BL_P) = -k(BL_P) .* Kzz(BL_P) .* ((h(BL_P+Nz) - h(BL_P))./dxe(BL_P));
flux_n_old(BL_P) = -k_old(BL_P) .* Kzz(BL_P) .* ((h_old(BL_P+1) - h_old(BL_P))./dzn(BL_P) + 1);
flux_e_old(BL_P) = -k_old(BL_P) .* Kzz(BL_P) .* ((h_old(BL_P+Nz) - h_old(BL_P))./dxe(BL_P));
F(BL_P) = psi(BL_P) - psi_old(BL_P) + ...
    (theta *dt *( (1./DX(BL_P)) .*(flux_e(BL_P))      +...
    (1./DZ(BL_P)) .* (flux_n(BL_P)))  +Q(BL_P)) + ...
    ((1-theta) *dt *( (1./DX(BL_P)) .*(flux_e_old(BL_P))      +...
    (1./DZ(BL_P)) .* (flux_n_old(BL_P)) )  +Q_old(BL_P)) ;
% 
%% BOT RIGHT
flux_n(BR_P) = -k(BR_P) .* Kzz(BR_P) .* ((h(BR_P+1) - h(BR_P))./dzn(BR_P) + 1);
flux_n_old(BR_P) = -k_old(BR_P) .* Kzz(BR_P) .* ((h_old(BR_P+1) - h_old(BR_P))./dzn(BR_P) + 1);
F(BR_P) = psi(BR_P) - psi_old(BR_P) + ...
    (theta *dt *( (1./DX(BR_P)) .*(Calc_RightBound(Z(BR_P),h(BR_P),Kxx(BR_P),CSG_HEIGHT) - flux_e(BR_P-Nz))      +...
    (1./DZ(BR_P)) .* (flux_n(BR_P)))  +Q(BR_P)) + ...
    ((1-theta) *dt *( (1./DX(BR_P)) .*(Calc_RightBound(Z(BR_P),h_old(BR_P),Kxx(BR_P),CSG_HEIGHT) - flux_e_old(BR_P-Nz))      +...
    (1./DZ(BR_P)) .* (flux_n_old(BR_P)) ) +Q_old(BR_P)) ;


 end