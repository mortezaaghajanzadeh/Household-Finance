%% Decomposition, using a sequence
clear all
main
% load('Results\main')  % can directly load solution; otherwise, run main.m
DepositProfit = solution0.Profit_deposit; save('Data\DepositProfit','DepositProfit')

a1_rec = Par.a1;
theta_rec = Par.theta;
ksi_rec = Par.ksi;


%%
% -  Reserve Requirement
Par.a1 = a1_rec;  Par.theta = 0.000; Par.ksi = ksi_rec;
Par = Grid_exp(Par); solution1A =solve_model_MP(Par,'YES DMP','YES LMP');
[firm1A, simu_con1A, simu_un1A] = SimulatePanel_short(Par, solution1A, shock0);

% - Capital Requirement
Par.a1 = a1_rec;  Par.theta = theta_rec; Par.ksi = 0.00;
Par = Grid_exp(Par); solution2A =solve_model_MP(Par,'YES DMP','YES LMP');
[firm2A, simu_con2A, simu_un2A] = SimulatePanel_short(Par, solution2A, shock0);

% - Financing Friction
Par.a1 = 0;  Par.theta = theta_rec; Par.ksi = ksi_rec;
Par = Grid_exp(Par); solution3A =solve_model_MP(Par,'YES DMP','YES LMP');
[firm3A, simu_con3A, simu_un3A] = SimulatePanel_short(Par, solution3A, shock0);

% - Deposit Market Power 
Par.a1 = a1_rec;  Par.theta = theta_rec; Par.ksi = ksi_rec;
Par = Grid_exp(Par); solution4A =solve_model_MP(Par,'NO DMP','YES LMP');
[firm4A, simu_con4A, simu_un4A] = SimulatePanel_short(Par, solution4A, shock0);

% Remove loan market power

Par.a1 = a1_rec;  Par.theta = theta_rec; Par.ksi = ksi_rec;
Par = Grid_exp(Par); solutionN = solve_model_noLMP(Par);

load('Data\shocks_6_500000_1_0_0'); 
[firmN, simu_conN, simu_unN] = SimulatePanel_short(Par, solutionN, shock0);

save("Results\T4_Determinants",'Par','simu_un0','simu_un1A','simu_un2A','simu_un3A','simu_un4A','simu_unN')
% -  Remove all


%% Print results
Par.a1 = 0;  Par.theta = 0.000; Par.ksi = 0;
Par = Grid_exp(Par); solution5A =solve_model_MP(Par,'YES DMP','YES LMP');
[firm5A, simu_con5A, simu_un5A] = SimulatePanel_short(Par, solution5A, shock0);
load("Results\T4_Determinants")

inc1 = (-simu_un1A.reg_B_agg + simu_un0.reg_B_agg)/simu_un0.reg_B_agg;
inc2 = (-simu_un2A.reg_B_agg + simu_un0.reg_B_agg)/simu_un0.reg_B_agg;
inc3 = (-simu_un3A.reg_B_agg + simu_un0.reg_B_agg)/simu_un0.reg_B_agg;
inc4 = (-simu_un4A.reg_B_agg + simu_un0.reg_B_agg)/simu_un0.reg_B_agg;
inc5 = (-simu_un5A.reg_B_agg + simu_un0.reg_B_agg)/simu_un0.reg_B_agg;
incN = (-simu_unN.reg_B_agg  + simu_un0.reg_B_agg)/simu_un0.reg_B_agg;

temp = [simu_un0.reg_B_agg, simu_un1A.reg_B_agg, simu_un2A.reg_B_agg, simu_un4A.reg_B_agg, simu_un5A.reg_B_agg, simu_unN.reg_B_agg;
    0, inc1, inc2, inc4, inc5, incN];

disp('Sensitivity  Change%')
disp(temp')
