%% Solve model under three cases: 1. baseline; 2. no capital regulation; and 3. deposit marekt has a "flat" profit
clear all
% main
load('Results\main')  % can directly load solution; otherwise, run main.m

DepositProfit = simu_un0.meanPdeposit;  save('Data\DepositProfit','DepositProfit')  

solution1 = solve_model_MP(Par,'NO DMP','YES LMP');
load('Data\shocks_6_500000_1_0_0'); 
[firm1, simu_con1, simu_un1, rec1] = SimulatePanel_short(Par, solution1, shock0);

Par.ksi = 0; Par = Grid_exp(Par); 
solution5 = solve_model_MP(Par,'YES DMP','YES LMP');
load('Data\shocks_6_500000_1_0_0'); 
[firm5, simu_con5, simu_un5, rec5] = SimulatePanel_short(Par, solution5, shock0);                                

save("Results\F5_PlotSolution", 'Par', 'simu_con0', 'simu_un0', 'simu_con1', 'simu_un1', 'simu_con5', 'simu_un5')                  

%% Plot
load("Results\F5_PlotSolution")
figure(5)

set(0, 'DefaultAxesFontSize', 12) 
set(0,'defaultTextInterpreter','latex');

subplot(2,2,1)
plot(Par.fGrid,smooth(simu_con0.L1/simu_con0.L1(6), 'lowess')); box off
title("Panel A: Bank lending",'FontSize', 12,'FontWeight', 'Normal');
xlabel('Federal funds rate')

subplot(2,2,2)
plot(Par.fGrid,smooth(simu_con5.L1/simu_con5.L1(6), 'lowess')); box off
title('Panel B: Unconstrained bank lending','FontSize', 12,'FontWeight', 'Normal');
xlabel('Federal funds rate')

subplot(2,2,3)
plot(Par.fGrid,smooth((simu_con0.profit+simu_con0.E0)/simu_con0.L1(6), 'lowess')); box off
title('Panel C: Bank capital','FontSize', 12,'FontWeight', 'Normal');
xlabel('Federal funds rate')

subplot(2,2,4)
plot(Par.fGrid,smooth(simu_con1.L1/simu_con1.L1(6), 'lowess')); box off
title('Panel D: Lending w/o deposit market power','FontSize', 12,'FontWeight', 'Normal');
xlabel('Federal funds rate')
ylim([0.8 1.45])


