%%
clear
close all
clc
useGPU = 0;
%%
nSimul = 100000;
%% Load training data
load("Allp.mat")
load('Allm.mat')
training.Allp = Allp;
training.Allm = Allm;

%% Estimation boundry
constrains.LB = [.8 1.01 0];
constrains.UB = [.99 5 .01];

%% Load Moments

tMoments = false(1,7); %%% just a test, include all

for v = [1 3 4 6] % 1-participation rate 3-conditional mean risky share 4-mean wealth 6-mean wealth at retirement
    tMoments(v) = true;
end
moment.tMoments = tMoments;
moment.targetm = targetm;
%% Estimation
x = point_estimate(training,moment,constrains);

m = simul([param.beta param.gamma param.phi],nSimul,useGPU,0.0);

disp("Estimated values")
disp(x)
disp("target moments:")
disp(moment.targetm)
disp("Estimated moments:")
disp(m)
