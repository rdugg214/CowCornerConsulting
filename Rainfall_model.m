clear
clc
%% Data collection
% Tipton
Tipton = readmatrix('Tipton.csv');
T_annual = Tipton(:,16).';
T_annual = [T_annual(1) NaN T_annual(2:end)];
T_time = Tipton(:,3).';
T_time = [T_time(1) 1969 T_time(2:end)];

% Bowenville
Bowenville = readmatrix('Bowenville.csv');
B_annual = Bowenville(:,16).';
B_time = Bowenville(:,3).';

% Cecil Plains
Cecil = readmatrix('Cecil_Plains.csv');
C_annual = Cecil(:,16);
C_time = Cecil(:,3).';

% Victory Downs
Victory = readmatrix('Victory_Downs.csv');
V_annual = Victory(:,16).';
V_time = Victory(:,3).';

% Plotting each location
% figure
% subplot(4,1,1)
% bar(T_time,T_annual)
% subplot(4,1,2)
% bar(B_time,B_annual)
% subplot(4,1,3)
% bar(C_time,C_annual)
% subplot(4,1,4)
% bar(V_time,V_annual)
%% Combining data into one annual vector

for i = 1:length(T_annual)
    if isnan(T_annual(i))
        year = T_time(i);
        k = find(B_time == year);
        T_annual(i) = B_annual(k);
    end
    if isnan(T_annual(i))
        k = find(C_time == year);
        T_annual(i) = C_annual(k);
    end
    if isnan(T_annual(i))
        k = find(V_time == year);
        T_annual(i) = V_annual(k);
    end
    
end

T_time = [1959:1967 T_time];
k1 = find(V_time == 1959);
k2 = find(V_time == 1967);

T_annual = [V_annual(k1:k2) T_annual];

%% Predicting future annual data
% Requried number of years
n = 200;
[prediction_data] = annual_predictor(T_annual,n);

%% Plotting predicted and previous data
figure
subplot(2,1,1)
bar(T_annual)
subplot(2,1,2)
bar(prediction_data)
save prediction_data prediction_data
%% Daily rainfall prediction
% Day
day = 1.9;
scale = 1; % normal rainfall 
[prediction, year] = rainfall(day,scale,prediction_data);


    



