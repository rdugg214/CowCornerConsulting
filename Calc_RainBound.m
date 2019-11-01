function fluxn_rain = Calc_RainBound(t,h,DX,rain,prediction_data)

locs = h<0;
fluxn_rain =zeros(size(h));
if rain == 1
    fluxn_rain(locs) = -(1 +  cos((2*pi*t)/365 ))*1.803 * 0.001 ;
elseif rain >=2
    tn = t+0.001;
    [ r, ~] = rainfall(1,tn,prediction_data);
    fluxn_rain(locs) = -r;
end


end