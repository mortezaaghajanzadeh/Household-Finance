%% Different Market Power
clear all
% main
load('Results\main')  % can directly load solution; otherwise, run main.m

Ngrid = [3,4,5,6,7,8,10];
for i = 1:length(Ngrid)
    Par.ND  = Ngrid(i);  Par.lD = 1.985 + 1.455 - log(Par.ND);  Par.nuB = 1.1509 - log(Par.NB) - Par.os;   Par.W0  = 0.2167*Par.ND/Par.NB; 
    Par = Grid_exp(Par); solution = solve_model_MP(Par,'YES DMP','YES LMP'); [~,~,simuNi, recNi] = SimulatePanel_short(Par, solution, shock0);
    tempND(i) = simuNi;
end
%%
main
for i = 1:length(Ngrid)
    Par.NB  = Ngrid(i);  Par.lD = 1.985 + 1.455 - log(Par.ND);  Par.nuB = 1.1509 - log(Par.NB) - Par.os;    Par.W0  = 0.2167*Par.ND/Par.NB; 
    Par = Grid_exp(Par); solution = solve_model_MP(Par,'YES DMP','YES LMP'); [~,~,simuNi, recNi] = SimulatePanel_short(Par, solution, shock0);
    tempNB(i) = simuNi;
end

%% Table
for i  = 1:length(Ngrid)
    rec(i,1) = tempNB(i).reg_D_agg;
    rec(i,2) = tempNB(i).reg_B_agg;
    rec(i,3) = tempNB(i).reg_rd;
    rec(i,4) = tempNB(i).reg_r;
    rec(i,5) = tempNB(i).SS_rate;
    rec(i,6) = tempNB(i).SS_quantity;

    rec(i,7) = tempND(i).reg_D_agg;
    rec(i,8) = tempND(i).reg_B_agg;
    rec(i,9) = tempND(i).reg_rd;
    rec(i,10) = tempND(i).reg_r;
    rec(i,11) = tempND(i).DSS_rate;
    rec(i,12) = tempND(i).DSS_quantity;
end
rec = [Ngrid',rec];

save("Results\F3_Concentration",'Ngrid','rec')

%% Plot
load("Results\F3_Concentration")

figure(3)

subplot(2,1,1)
set(0, 'DefaultAxesFontSize', 11) 
yyaxis left
temp = smooth(rec(:,12)); plot(1./Ngrid,temp); % change relative to the benchmark
ylabel("Sensitivity to federal funds rate"+newline+"   ",'interpreter','latex')
yyaxis right
temp = smooth(rec(:,13)); plot(1./Ngrid,temp);
xlabel("Bank concentration"+newline+"   ",'interpreter','latex')
xlim([0.1, 0.33])
legend("Deposit rate","Deposit quantity",'Orientation','horizontal','location','southoutside','interpreter','latex'); legend boxoff;

%
subplot(2,1,2)
set(0, 'DefaultAxesFontSize', 11) 
yyaxis left
temp = smooth(rec(:,6)); plot(1./Ngrid,temp);
ylabel("Sensitivity to federal funds rate"+newline+"   ",'interpreter','latex')
yyaxis right
temp = smooth(rec(:,7)); plot(1./Ngrid,temp);
xlabel("Bank concentration"+newline+"   ",'interpreter','latex')
xlim([0.1, 0.33])
legend("Loan rate","Loan quantity",'Orientation','horizontal','location','southoutside','interpreter','latex'); legend boxoff;

saveas(gcf,'replication_figure4.eps','epsc')

