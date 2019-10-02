function fluxn_rain = Calc_RainBound(t,h)

if h<0
    fluxn_rain = 1.793 + 0.994 * cos((2*pi*t)/365);
else
    fluxn_rain = 0;
end
end