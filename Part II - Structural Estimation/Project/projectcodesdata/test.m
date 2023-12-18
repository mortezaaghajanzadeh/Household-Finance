
clear
close all
clc



useGPU = 0;
%%
nSimul = 100000;
%%% foundamentals [to be estimated, random numbers here: as an example]
param.beta = 0.92;
param.gamma = 2;
param.phi = .005;
m = simul([param.beta param.gamma param.phi],nSimul,useGPU,0.0);

%%
m2 = simul_counter([param.beta param.gamma param.phi],nSimul,useGPU,0.0);

%%

%%

clear
close all
clc



useGPU = 0;

%%
load("Allp.mat")
load('Allm.mat')
training.Allp = Allp;
training.Allm = Allm;
%%

nSimul = 100000;
constrains.LB = [.8 1.01 0];
constrains.UB = [.99 5 .01];

%%
clc


tMoments = false(1,7); %%% just a test, include all

for v = [1 3 4 6]
    tMoments(v) = true;
end

targetm = readtable('data_moments.txt'); targetm = targetm{:,:};

moment.tMoments = tMoments;
moment.targetm = targetm;

%% 
x = point_estimate(training,moment,constrains)