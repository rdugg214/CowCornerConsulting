% parameters 
clear, close all, clc
dtmax = 3720; % max time step
endtime = 20*365; 
t_on_CSG = 20*365; %time on of CSG 
t_on_PUMP = 20*365; %time on of PUMP
Pr = 0.25; % Pumping rate 
dx =50; dz =10; %distance between nodes
nx = 21; ny = 11;  %number of nodes
geometric = [dx dz]; %uniform progression
% geometric = [1 nx ny]; %geometric progression
simple = 1; %1 if river, evap and complex rain on.
SAVEVID = 1; %Set 1 if want a video to save, 0 otherwise

%default simulation parameters
DEF = v2struct(dtmax, endtime, t_on_CSG,t_on_PUMP,geometric, simple,Pr,SAVEVID);

%%
SET = DEF;
SET.simple =1 ;
h = SimpleSolution(SET)
