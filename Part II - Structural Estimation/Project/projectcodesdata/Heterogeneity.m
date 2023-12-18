%%
clear
close all
clc
useGPU = 0;
%%
disp("The estimated simulation")
nSimul = 100000;
%%% foundamentals [to be estimated, random numbers here: as an example]
param.beta = 0.9519;
param.gamma = 2.0143;
param.phi = .0024;
param.capital_tax = 0.0;
m = simul([param.beta param.gamma param.phi param.capital_tax],nSimul,useGPU);
%% Type I
disp("The simulation for type I")
nSimul_1 = nSimul/2;
param_1.beta = param.beta + 0.01;
param_1.gamma = param.gamma - 0.5;
param_1.phi = param.phi * 0;
param_1.capital_tax = 0.0;
m_1 = simul([param_1.beta param_1.gamma param_1.phi],nSimul_1,useGPU);
%% Type II
disp("The simulation for type II")
nSimul_2 = nSimul/2;
param_2.beta = param.beta - 0.01;
param_3.gamma = param.gamma + 0.5;
param_2.phi = param.phi * 2;
param_2.capital_tax = 0.0;
m_2 = simul([param_2.beta param_2.gamma param_2.phi],nSimul_2,useGPU);

%% Load training data
load("Allp.mat")
load('Allm.mat')
training.Allp = Allp;
training.Allm = Allm;

%% Estimation boundry
constrains.LB = [.8 1.01 0];
constrains.UB = [.99 5 .01];

%% Moment selection

tMoments = false(1,7); %%% just a test, include all

for v = [1 3 4 6] % 1-participation rate 3-conditional mean risky share 4-mean wealth 6-mean wealth at retirement
    tMoments(v) = true;
end
moment.tMoments = tMoments;
moment.targetm = (m_1 * nSimul_1 + m_2 * nSimul_2) / (nSimul_1 + nSimul_2);

disp("The estimation by hetrogenous moments")
x = point_estimate(training,moment,constrains);
