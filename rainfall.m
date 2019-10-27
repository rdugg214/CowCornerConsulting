function [prediction, k] = rainfall(scale,day,predictor_data)
k = floor(day/365) + 1;
D = predictor_data(k)/365;

d = mod(day,365);
A = 0.9944;
deviation = 1.1714;
prediction = scale*deviation * A * cos(2*pi*d/365) + D;

if prediction < 0
   prediction = 0;
end

prediction = prediction/1000;
end

