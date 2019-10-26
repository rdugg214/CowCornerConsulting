% parameters 
clear, close all, clc
dtmax = 360; % max time step
endtime = 100*365; 
t_on_CSG = 5*365; %time on of CSG 
t_on_PUMP = 5*365; %time on of PUMP
Pr = 0.25; % Pumping rate 
dx =50; dz = 5; %distance between nodes
nx = 21; ny = 11;  %number of nodes
geometric = [dx dz]; %uniform progression
% geometric = [1 nx ny]; %geometric progression
simple = 1; %river and evap off

DEF = v2struct(dtmax, endtime, t_on_CSG,t_on_PUMP,geometric, simple,Pr);
% %% Run at uniform only rain and pumping
% SET = DEF;
% SET.endtime = 2*365;
% SimpleSolution(SET)
%%
SET = DEF;
SET.simple = 1;
SET.t_on_CSG = 2*365;
SimpleSolution(SET)
