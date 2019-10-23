function fluxn_rain = Calc_RainBound(t,h,DX)


fluxn_rain =zeros(size(h));
locs = h<0;
    fluxn_rain(locs) = -(1 +  cos((2*pi*t)/365 ))*1.803 * 0.001 ;

% fluxn_rain(locs) = -0.003;
% fint =   @(t) (t.*(-1.803e-3)-(sin(t.*pi.*(2.0./3.65e+2)).*3.290475e-1)./pi);
% f_av_yr = 0.6544;
% fav_day = 0.0018;

end