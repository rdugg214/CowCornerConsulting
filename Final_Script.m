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
Pr = 0.25; % Pumping rate
dx =50; dz =2; %distance between nodes
nx = 11; nz = 31;  %number of nodes
geometric = [dx dz]; %uniform progression
% geometric = [1 nx nz]; %geometric progression
SAVEVID = 1; %Set 1 if want a video to save, 0 otherwise
t_on_CSG = 2000*365; %time on of CSG 
t_on_PUMP = 2000*365; %time on of PUMP
DELCSG = 5000;
simple = 1; %1 if river, evap and complex rain on.
rain = 2; %0 - no rain, 1 - simple, 2 - complex rain, 3-complex drought rain
h_init = [];
t_init = 0;
NAME = 'DEF';
% r = load('RES.mat');
% h_init = r.RES.h_final;
% t_init = r.RES.t_final;
DEF = v2struct(dtmax, endtime, t_on_CSG,t_on_PUMP,geometric, simple,Pr,SAVEVID,DELCSG,h_init,t_init,rain,NAME);
run_sim_Rain_Test = 0;
run_sim_full_run = 0;
run_sim_normal_csg0 = 0;
run_sim_normal_csg1 = 0;
run_sim_normal_csg2 = 0;
run_sim_normal_csg3 = 0;
run_sim_drought_csg0 = 0;
run_sim_drought_csg1 = 0;
run_sim_drought_csg2 = 0;
run_sim_drought_csg3 = 0;



%% Compare with Analytical Rain and Pumping
if run_sim_Rain_Test
SET = DEF;
SET.simple =1;
SET.rain = 1;
SET.t_on_CSG = Inf;
SET.t_on_PUMP = 0;
SET.endtime = 10*365;
SET.SAVEVID = 1;
SET.NAME = 'RES_Rain_Test';
% Run Simulation
RES_Rain_Test = Final_Simulation(SET);
save('RES_Rain_Test.mat','RES_Rain_Test')
else
load('RES_Rain_Test.mat')
end

%% Running simulation normal rainfall: Running up to 20 years
%default simulation parameters

if run_sim_full_run
SET = DEF;
SET.simple =0;
SET.endtime = 21*365;
SET.SAVEVID = 1;
RES_full_run = Final_Simulation(SET);
SET.NAME = 'RES_full_run';
RES = RES_rain_nocsg;
save('RES.mat','RES')
else
load('RES.mat')
end


%% Running simulation normal rainfall no CSG
if run_sim_normal_csg0
SET = DEF;
SET.h_init = 1;
SET.t_init = 20*365;
SET.endtime = 50*365;
SET.t_on_CSG = Inf;
SET.t_on_PUMP = 0*365;
SET.DELCSG = Inf; 
SET.simple =0;
SET.NAME = 'RES_rain_csg0';
%default simulation parameters
RES_rain_csg0 = Final_Simulation(SET);
save('RES_rain_csg0.mat','RES_rain_csg0')
else
load('RES_rain_csg0.mat') 
end

%% Running simulation normal rainfall at location 1
if run_sim_normal_csg1
SET = DEF;
r = load('RES.mat');
SET.h_init = 1;
SET.t_init = 20*365;
SET.endtime = 50*365;
SET.t_on_CSG = 0*365;
SET.t_on_PUMP = 0*365;
SET.DELCSG = 1000; 
SET.simple =0;
SET.NAME = 'RES_rain_csg1';
%default simulation parameters
RES_rain_csg1 = Final_Simulation(SET);
save('RES_rain_csg1.mat','RES_rain_csg1')
else
load('RES_rain_csg1.mat') 
end
%% Running simulation normal rainfall at location 2
if run_sim_normal_csg2
SET = DEF;
r = load('RES.mat');
SET.h_init = 1;
SET.t_init = 20*365;
SET.endtime = 50*365;
SET.t_on_CSG = 0*365;
SET.t_on_PUMP = 0*365;
SET.DELCSG = 2500; 
SET.simple =0;
SET.NAME = 'RES_rain_csg2';
%default simulation parameters
RES_rain_csg2 = Final_Simulation(SET);
save('RES_rain_csg2.mat','RES_rain_csg2')
else
load('RES_rain_csg2.mat') 
end
%% Running simulation normal rainfall at location 3
if run_sim_normal_csg3
SET = DEF;
r = load('RES.mat');
SET.h_init = 1;
SET.t_init = 20*365;
SET.endtime = 50*365;
SET.t_on_CSG = 0*365;
SET.t_on_PUMP = 0*365;
SET.DELCSG = 5000; 
SET.simple =0;
SET.NAME = 'RES_rain_csg3';
%default simulation parameters
RES_rain_csg3 = Final_Simulation(SET);
save('RES_rain_csg3.mat','RES_rain_csg3')
else
load('RES_rain_csg3.mat') 
end

%% Running simulation drought rainfall no CSG
if run_sim_drought_csg0
close all
SET = DEF;
r = load('RES.mat');
SET.h_init =1;
SET.t_init = 20*365;
SET.endtime = 50*365;
SET.t_on_CSG = Inf;
SET.rain = 3;
SET.t_on_PUMP = 0*365;
SET.DELCSG = Inf; 
SET.simple =0;
SET.NAME = 'RES_drought_csg0';
%default simulation parameters
RES_drought_csg0 = Final_Simulation(SET);
save('RES_drought_csg0.mat','RES_drought_csg0')
else
load('RES_drought_csg0.mat') 
end
% 
%% Running simulation drought rainfall at location 1
if run_sim_drought_csg1
SET = DEF;
r = load('RES.mat');
SET.h_init = 1;
SET.t_init = 20*365;
SET.endtime = 50*365;
SET.t_on_CSG = 0*365;
SET.rain = 3;
SET.t_on_PUMP = 0*365;
SET.DELCSG = 1000; 
SET.simple =0;
SET.NAME = 'RES_drought_csg1';
%default simulation parameters
RES_drought_csg1 = Final_Simulation(SET);
save('RES_drought_csg1.mat','RES_drought_csg1')
else
load('RES_drought_csg1.mat') 
end

%% Running simulation drought rainfall at location 2
if run_sim_drought_csg2
SET = DEF;
r = load('RES.mat');
SET.h_init = 1;
SET.t_init = 20*365;
SET.endtime = 50*365;
SET.t_on_CSG = 0*365;
SET.rain = 3;
SET.t_on_PUMP = 0*365;
SET.DELCSG = 2500; 
SET.simple =0;
SET.NAME = 'RES_drought_csg2';
%default simulation parameters
RES_drought_csg2 = Final_Simulation(SET);
save('RES_drought_csg2.mat','RES_drought_csg2')
else
load('RES_drought_csg2.mat') 
end

%% Running simulation drought rainfall at location 3
if run_sim_drought_csg3
SET = DEF;
r = load('RES.mat');
SET.h_init = 1;
SET.t_init = 20*365;
SET.endtime = 50*365;
SET.t_on_CSG = 0*365;
SET.t_on_PUMP = 0*365;
SET.rain = 3;
SET.DELCSG = 5000; 
SET.simple =0;
SET.NAME = 'RES_drought_csg3';
%default simulation parameters
RES_drought_csg3 = Final_Simulation(SET);
save('RES_drought_csg3.mat','RES_drought_csg3')
else
load('RES_drought_csg3.mat') 
end


