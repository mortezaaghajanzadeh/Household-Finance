%%

clear
close all
clc

useGPU = 0;

%% Parameters
load("estimated_point")
param.beta = x(1);
param.gamma = x(2);
param.phi = x(3);
nSimul = 100000;
%% baseline
m = simul([param.beta param.gamma param.phi],nSimul,useGPU,0.0);
%% adding tax
m_tax = simul([param.beta param.gamma param.phi],nSimul,useGPU,0.3);
%% Add tax with zero participation cost
m_0 = simul([param.beta param.gamma 0],nSimul,useGPU,0.3);

%% Add tax with more risk aversion
m_risk = simul([param.beta param.gamma * 2 param.phi],nSimul,useGPU,0.3);

%% results
results = [m; m_tax; m_0; m_risk];
T = array2table(results');

T.Properties.VariableNames(1:4) = {'$\tau =0$','$\tau =30\%$','$\tau =30\% \& \phi = 0$','$\tau =30\% \& \gamma = 2\gamma$'};
T.Properties.RowNames(1:7) = ["$E[\pi_{it}>0]$", "$E[\pi_{it}]$","$E[\pi_{it}|\pi_{it}>0]$", '$E[W_{it}]$',"$\sigma[W_{it}]$",'$E[W_{iT}]$',"$\sigma[W_{iT}]$"];
table2latex(T,'./counterfactual_moments.tex')