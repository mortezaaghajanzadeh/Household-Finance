%%

clear
close all
clc

useGPU = 0;

%% Parameters
param.beta = 0.9519;
param.gamma = 2.0143;
param.phi = .0024;

%% baseline
nSimul = 100000;
m = simul([param.beta param.gamma param.phi],nSimul,useGPU,0.0);

%% adding tax
m_tax = simul([param.beta param.gamma param.phi],nSimul,useGPU,0.3);

%% Add tax with zero participation cost
m_0 = simul([param.beta param.gamma 0],nSimul,useGPU,0.3);

%% Add tax with more risk aversion
m_risk = simul([param.beta param.gamma * 2 param.phi],nSimul,useGPU,0.3);
