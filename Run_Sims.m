% parameters 
clear, close all, clc
dtmax = 100; % max time step
endtime = 20*365; 
t_on_CSG = 30*365; %time on of CSG 
t_on_PUMP = 30*365; %time on of PUMP
Pr = 0.5; % Pumping rate 
DELCSG = 5000;
dx =50; dz =10; %distance between nodes
nx = 21; ny = 11;  %number of nodes
geometric = [dx dz]; %uniform progression
% geometric = [1 nx ny]; %geometric progression
simple = 1; %1 if river, evap and complex rain on.
SAVEVID = 1; %Set 1 if want a video to save, 0 otherwise
h_init = [];
t_init = 0;

% r = load('RES.mat');
% h_init = r.RES.h_final;
% t_init = r.RES.t_final;

%default simulation parameters
DEF = v2struct(dtmax, endtime, t_on_CSG,t_on_PUMP,geometric, simple,Pr,SAVEVID,DELCSG,h_init,t_init);

%%
SET = DEF;
SET.simple =0;
RES = SimpleSolution(SET)
%%
h = RES.h_final;
%% save
% save('RES.mat','RES')