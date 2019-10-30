function qw_river = Calc_RiverBound(z,hp,simple)

qw_river = 0;
if ~simple
HR = 85;
KR = 0.1;
Rt = 30;
Rb = 80;
Rd = 20;
if hp + z <= HR
       qw_river = KR*(HR-(hp+z))/Rt;
elseif hp + z < Rb
       qw_river = -KR*(HR-Rb)/Rt;
elseif HR < (hp + z) && (hp + z) <=100
       qw_river = KR*(HR-(hp+z))/Rt;
else
qw_river = 0;
end
end
% disp('river error')
% end

end