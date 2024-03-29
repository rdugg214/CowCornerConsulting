%% Cow Corner Consulting Final Script
% Running the script outputs all results found in the report.
% The script does take a long duration to run, prepare for long
% wait after running
%
% clear, close all, clc
% 
%% Defining Consistant Parameters for all simulations

clear, close all, clc
dtmax = 180; % max time step
endtime = 21*365;
Pr = 0.5; % Pumping rate
dx =50; dz =2; %distance between nodes
nx = 11; nz = 31;  %number of nodes
% geometric = [dx dz]; %uniform progression
geometric = [1 nx nz]; %geometric progression
SAVEVID = 0; %Set 1 if want a video to save, 0 otherwise
t_on_CSG = 2000*365; %time on of CSG 
t_on_PUMP = 2000*365; %time on of PUMP
DELCSG = 5000;
simple = 1; %1 if river, evap and complex rain on.
rain = 2; %0 - no rain, 1 - simple, 2 - complex rain, 3-complex drought rain
h_init = [];
t_init = 0;
% r = load('RES.mat');
% h_init = r.RES.h_final;
% t_init = r.RES.t_final;
DEF = v2struct(dtmax, endtime, t_on_CSG,t_on_PUMP,geometric, simple,Pr,SAVEVID,DELCSG,h_init,t_init,rain);
%% Compare with Analytical Rain and Pumping
SET = DEF;
% r = load('RES.mat');
% SET.h_init = [];
% SET.t_init = r.RES.t_final;
SET.geometric = [50 5];
SET.simple =1;
SET.t_on_PUMP = 40*365;
SET.t_on_PUMP = 40*365;
SET.endtime = 40*365;
SET.rain = 1;
SET.SAVEVID = 0;
% Run Simulation
RES_Rain_Test = Final_Simulation(SET);
%%
% 
%% Running simulation normal  rainfall without CSG
% 
% %default simulation parameters
% SET = DEF;
% SET.simple =0;
% % Run Simulation
% RES_rain_nocsg = Final_Simulation(SET);
% RES = RES_rain_nocsg;
% save('RES.mat','RES')
% Running simulation normal rainfall with no CSG
close all
SET = DEF;
r = load('RES.mat');
% SET.h_init = r.RES.h_final;
SET.t_init = r.RES.t_final;
SET.endtime = 50*365;
SET.t_on_CSG = 50*365;
SET.t_on_PUMP = 50*365;
SET.Pr= 0.25;
SET.DELCSG = 1000; 
% SET.dtmax = 20;
SET.simple =0;
%default simulation parameters
RES_rain_csg0 = Final_Simulation(SET);
% %% Running simulation normal rainfall with CSG at distance 1
% close all
% SET = DEF;
% r = load('RES.mat');
% % SET.h_init = r.RES.h_final;
% SET.t_init = r.RES.t_final;
% SET.endtime = 50*365;
% SET.t_on_CSG = 0*365;
% SET.t_on_PUMP = 50*365;
% SET.Pr= 0.25;
% SET.DELCSG = 1000; 
% SET.dtmax = 20;
% SET.simple =0;
% %default simulation parameters
% 
% % Run Simulation
% RES_rain_csg1 = Final_Simulation(SET);

% %% Running simulation normal rainfall with CSG at distance 2
% close all
% SET = DEF;
% r = load('RES.mat');
% % SET.h_init = r.RES.h_final;
% SET.t_init = r.RES.t_final;
% SET.endtime = 50*365;
% SET.t_on_CSG = 0*365;
% SET.t_on_PUMP =  50*365;
% SET.Pr= 0.25;
% SET.DELCSG = 5000; 
% SET.dtmax = 20;
% SET.simple =0;
% %default simulation parameters
% 
% 
% % Run Simulation
% RES_rain_csg2 = Final_Simulation(SET);

%% Running simulation normal rainfall with CSG at distance 3
close all
SET = DEF;
r = load('RES.mat');
% SET.h_init = r.RES.h_final;
SET.t_init = r.RES.t_final;
SET.endtime = 50*365;
SET.t_on_CSG = 0*365;
SET.t_on_PUMP = 50*365;
SET.Pr= 0.25;
SET.DELCSG = 10000; 
SET.dtmax = 20;
SET.simple =0;
%default simulation parameters


% Run Simulation
RES_rain_csg3 = Final_Simulation(SET);

% %% Running simulation drought rainfall with CSG at distance 1
% t_on_CSG = 30*365; %time on of CSG 
% t_on_PUMP = 30*365; %time on of PUMP
% DELCSG = 5000;
% simple = 1; %1 if river, evap and complex rain on.
% h_init = [];
% t_init = 0;
% 
% % r = load('RES.mat');
% % h_init = r.RES.h_final;
% % t_init = r.RES.t_final;
% 
% %default simulation parameters
% DEF = v2struct(dtmax, endtime, t_on_CSG,t_on_PUMP,geometric, simple,Pr,SAVEVID,DELCSG,h_init,t_init);
% SET = DEF;
% SET.simple =0;
% % Run Simulation
% RES_drought_csg1 = Final_Simulation(SET);

%% Running simulation normal rainfall with CSG at distance 2

%% Running simulation drought rainfall with CSG at distance 2

%% Running simulation normal rainfall with CSG at distance 3

%% Running simulation drought rainfall with CSG at distance 3

%% Plotting results for Analysis
% %%%%%%%%%%%%%%%%%%%%%%% Example Data
% norm_rain_no_csg_H = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10; 100, 79, 80, 84, 84, 80, 81, 82, 83, 89, 95];
% norm_rain_csg_1_H = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10; 100, 79, 80, 84, 84, 80, 81, 60, 68, 69, 69];
% drought_no_csg_H = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10; 100, 79, 80, 84, 84, 80, 81, 70, 70, 75, 75];
% drought_csg_1_H = [0, 0.1, 2, 3, 6, 6.5, 6.9, 7, 8, 9; 100, 79, 80, 84, 84, 80, 81, 70, 40, 20];
% 
% figure();
% hold on;
% plot(norm_rain_no_csg_H(1, :), norm_rain_no_csg_H(2, :), 'r-');
% plot(norm_rain_csg_1_H(1, :), norm_rain_csg_1_H(2, :), 'r--');
% plot(drought_no_csg_H(1, :), drought_no_csg_H(2, :), 'b-');
% plot(drought_csg_1_H(1, :), drought_csg_1_H(2, :), 'b--');
% legend("Normal Rainfall - No CSG", "Normal Rainfall - CSG 1", "Drought - No CSG", "Drought - CSG 1")
% title("Water table height over time");
% xlabel("Time (Days)");
% ylabel("Water table height (m)")
% %%%%%%%%%%%%%%%%%%%%%%%


% figure();
% hold on;
% plot(RES_rain_nocsg.t_vector, RES_rain_nocsg.Ballarr(1, :), 'r-');
% plot(RES_rain_csg1.t_vector, RES_rain_csg1.Ballarr(1, :), 'r--');
% plot(RES_drought_nocsg.t_vector, RES_drought_nocsg.Ballarr(1, :), 'b-');
% plot(RES_drought_csg1.t_vector, RES_drought_csg1.Ballarr(1, :), 'b--');
% legend("Normal Rainfall - No CSG", "Normal Rainfall - CSG 1", "Drought - No CSG", "Drought - CSG 1")
% title("Flux From River");
% xlabel("Time (Days)");
% ylabel("Flux")
clear
Lx = 500;
Lz = 100;

load('RES_drought_csg0.mat')
load('RES_drought_csg1.mat')
load('RES_drought_csg2.mat')
load('RES_drought_csg3.mat')
load('RES_rain_csg0.mat')
load('RES_rain_csg1.mat')
load('RES_rain_csg2.mat')
load('RES_rain_csg3.mat')
% RES_drought_csg3.Harr = RES_drought_csg3.Harr./(Lx*Lz).*length(RES_drought_csg3.Harr)
% save('RES_drought_csg3.mat','RES_drought_csg3')
%%
load('RES.mat')
close all
figure();
hold on;
% plot(RES_rain_nocsg.t_vector, RES_rain_nocsg.Harr, 'r-');

plot(RES_rain_csg1.t_vector+RES.t_final, RES_rain_csg1.Harr, 'r','LineWidth',1);
plot(RES_rain_csg2.t_vector+RES.t_final, RES_rain_csg2.Harr, 'b','LineWidth',3);
plot(RES_rain_csg3.t_vector+RES.t_final, RES_rain_csg3.Harr, 'g','LineWidth',3);
plot(RES_rain_csg0.t_vector+RES.t_final, RES_rain_csg0.Harr, 'k','LineWidth',4);
legend("CSG 1km", "CSG 5km", "CSG 10km","No CSG")
title("Water table height over time for Normal Rainflal");
xlabel("Time (Days)");
ylabel("Water table height (m)")

figure();
hold on;
% plot(RES_rain_nocsg.t_vector, RES_rain_nocsg.Harr, 'r-');

plot(RES_drought_csg1.t_vector+RES.t_final, RES_drought_csg1.Harr, 'r','LineWidth',1);
plot(RES_drought_csg2.t_vector+RES.t_final, RES_drought_csg2.Harr, 'b','LineWidth',3);
plot(RES_drought_csg3.t_vector+RES.t_final, RES_drought_csg3.Harr, 'g','LineWidth',3);
plot(RES_drought_csg0.t_vector+RES.t_final, RES_drought_csg0.Harr, 'k','LineWidth',4);
legend("CSG 1km", "CSG 5km", "CSG 10km","No CSG")
title("Water table height over time for Drought Rainflal");
xlabel("Time (Days)");
ylabel("Water table height (m)")

figure();
hold on;
% plot(RES_rain_nocsg.t_vector, RES_rain_nocsg.Harr, 'r-');

plot(RES_drought_csg1.t_vector+RES.t_final, -RES_drought_csg1.Ballarr(1,:), 'r','LineWidth',2);
% plot(RES_drought_csg1.t_vector+RES.t_final, -RES_drought_csg1.Ballarr(3,:), 'r-','LineWidth',2);
plot(RES_drought_csg2.t_vector+RES.t_final, -RES_drought_csg2.Ballarr(1,:), 'b','LineWidth',3);
plot(RES_drought_csg3.t_vector+RES.t_final, -RES_drought_csg3.Ballarr(1,:), 'g','LineWidth',3);
plot(RES_drought_csg0.t_vector+RES.t_final, -RES_drought_csg0.Ballarr(1,:), 'k','LineWidth',4);
legend("CSG 1km", "CSG 5km", "CSG 10km","No CSG")
title("Water flow INTO river for drought");
xlabel("Time (Days)");
ylabel("Flux (m/d)")
