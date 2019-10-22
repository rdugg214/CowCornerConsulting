% parameters 
dtmax = 360; % max time step
endtime = 100*365; 
t_on_CSG = 365*10; %time on of CSG 
t_on_PUMP = 365*10; %time on of PUMP

dx =50; dz = 2; %distance between nodes
nx = 21; ny = 11;  %number of nodes
geometric = [dx dz]; %uniform progression
% geometric = [1 nx ny]; %geometric progression
SimpleSolution(dtmax, endtime, t_on_CSG,t_on_PUMP,geometric)
