function qw_river = Calc_RiverBound(z,hp,simple,prediction_data,t)


if ~simple
    HR = 85;
    KR = 0.1;
    Rt = 30;
    Rb = 80;
    Rd = 20;
    if  hp + z <= HR && prediction_data(floor(t/365) + 1)<100
        qw_river = 0;
    elseif  hp + z < Rb&& prediction_data(floor(t/365) + 1)<100
        qw_river = 0;
    elseif hp + z <= HR
        qw_river = KR*(HR-(hp+z))/Rt;
    elseif hp + z < Rb
        qw_river = -KR*(HR-Rb)/Rt;
    elseif HR < (hp + z) && (hp + z) <=100
        qw_river = KR*(HR-(hp+z))/Rt;
    else
        qw_river = 0;
    end
else
    qw_river = 0;
end
% disp('river error')
% end

end