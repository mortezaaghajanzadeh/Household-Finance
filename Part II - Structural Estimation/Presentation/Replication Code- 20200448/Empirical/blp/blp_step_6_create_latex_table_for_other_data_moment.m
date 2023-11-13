clear
clc
filepath = '.\intermediate_output\';
texfolder=strcat('../','output\Tables');

i = 1; 

%fprintf(FID, '\\documentclass[12pt]{article}\\begin{document}');

% these moments are the same across sample
[~,~,payout_vw_1] = xlsread(strcat(filepath,'moment_payout_vw_1','.csv'));
[~,~,VAR] = xlsread(strcat(filepath,'VAR','.xlsx'));
[~,~,quality_dep] = xlsread(strcat(filepath,'demand_dep_',num2str(i),'.xls'));
[~,~,quality_loan] = xlsread(strcat(filepath,'demand_loan_',num2str(i),'.xls'));


for i = 1:5
 FID = fopen(strcat(texfolder,'\moment_',num2str(i),'.tex'), 'w');
% these moments change across sample
[~,~,demand_dep] = xlsread(strcat(filepath,'demand_dep_',num2str(i),'.xls'));
[~,~,demand_loan] = xlsread(strcat(filepath,'demand_loan_',num2str(i),'.xls'));
[~,~,spread_FFR_sensitivity] = xlsread(strcat(filepath,'moment_spread_FFR_sensitivity_',num2str(i),'.csv'));
[~,~,balance_sheet] = xlsread(strcat(filepath,'moment_bank_vw_',num2str(i),'.csv'));
[~,~,HHI] = xlsread(strcat(filepath,'HHI_',num2str(i),'.csv'));
if i<=3; [~,~,payout_vw_1] = xlsread(strcat(filepath,'moment_payout_vw_',num2str(i),'.csv'));end
    

fprintf(FID, '\\begin{tabular*}{\\hsize}{@{\\hskip\\tabcolsep\\extracolsep\\fill}l*{2}{c}}');
fprintf(FID, '\\hline ');
fprintf(FID, '\\hline ');
fprintf(FID, '& {Actual Moment}  &    \\\\ [1ex] \\hline  ');

%fprintf(FID, strcat('Sample',num2str(i)));
%fprintf(FID, '&');
%fprintf(FID, '\\\\\n');
fprintf(FID, 'Dividend yield&');
fprintf(FID, num2str(payout_vw_1{ismember(payout_vw_1(:,2),{'total_payout_ratio_w'}),1},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
if i<=3
    fprintf(FID, 'MB&');
    fprintf(FID, num2str(payout_vw_1{ismember(payout_vw_1(:,2),{'MB_w'}),1},'%.3f'));
    fprintf(FID, '&');
    fprintf(FID, '\\\\\n');
end
fprintf(FID, 'Leverage&');
fprintf(FID, num2str(balance_sheet{ismember(balance_sheet(:,2),{'leverage'}),1},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
% fprintf(FID, 'Borrowing/Assets&');
% fprintf(FID, num2str(balance_sheet{ismember(balance_sheet(:,2),{'borrowing_assets'}),1},'%.3f'));
% fprintf(FID, '&');
% fprintf(FID, '\\\\\n');
% fprintf(FID, 'Borrowing/Assets (sd)&');
% fprintf(FID, num2str(balance_sheet{ismember(balance_sheet(:,2),{'borrowing_assets_sd'}),1},'%.3f'));
% fprintf(FID, '&');
% fprintf(FID, '\\\\\n');
% fprintf(FID, 'Borrowing/Assets (demeaned sd)&');
% fprintf(FID, num2str(balance_sheet{ismember(balance_sheet(:,2),{'borrowing_assets_demean_sd'}),1},'%.3f'));
% fprintf(FID, '&');
% fprintf(FID, '\\\\\n');
fprintf(FID, 'Borrowing/Deposits&');
fprintf(FID, num2str(balance_sheet{ismember(balance_sheet(:,2),{'borrowing_dep_total'}),1},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, 'Borrowing/Deposits (sd)&');
fprintf(FID, num2str(balance_sheet{ismember(balance_sheet(:,2),{'borrowing_dep_total_sd'}),1},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, 'Borrowing/Deposits (demeaned sd)&');
fprintf(FID, num2str(balance_sheet{ismember(balance_sheet(:,2),{'borrowing_dep_total_demean_sd'}),1},'%.3f'));
% fprintf(FID, '&');
% fprintf(FID, '\\\\\n');
% fprintf(FID, 'Loan to Deposits Ratio&');
% fprintf(FID, num2str(balance_sheet{ismember(balance_sheet(:,2),{'loan_to_dep'}),1},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, 'Deposits to Assets Ratio&');
fprintf(FID, num2str(balance_sheet{ismember(balance_sheet(:,2),{'deposits_assets'}),1},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, 'Net noninterest expenses&');
fprintf(FID, num2str(-balance_sheet{ismember(balance_sheet(:,2),{'nonintinc_assets'}),1}+balance_sheet{ismember(balance_sheet(:,2),{'nonintexp_assets'}),1},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, 'Deposit Spreads&');
fprintf(FID, num2str(balance_sheet{ismember(balance_sheet(:,2),{'s_dep'}),1},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, 'Loan Spreads&');
fprintf(FID, num2str(balance_sheet{ismember(balance_sheet(:,2),{'s_loan'}),1}-balance_sheet{ismember(balance_sheet(:,2),{'netloanchargeoffs_loans'}),1},'%.3f')); %adjust expected default
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, 'Deposit Spreads - FFR Sensitivity&');
fprintf(FID, num2str(spread_FFR_sensitivity{ismember(spread_FFR_sensitivity(:,2),{'s_dep_FFR'}),1},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, 'Loan Spreads - FFR Sensitivity&');
fprintf(FID, num2str(spread_FFR_sensitivity{ismember(spread_FFR_sensitivity(:,2),{'s_loan_FFR'}),1},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, 'Total Credit - FFR Sensitivity&');
fprintf(FID, num2str(VAR{ismember(VAR(:,1),{'TotalCredit-FFR 3-year sensitivity: '}),2},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, 'Total Credit - FFR Sensitivity (se)&');
fprintf(FID, num2str(VAR{ismember(VAR(:,1),{'TotalCredit-FFR 3-year sensitivity: '}),3},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, 'Bank Loan - FFR Sensitivity&');
fprintf(FID, num2str(VAR{ismember(VAR(:,1),{'BankLoan-FFR 3-year sensitivity: '}),2},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, 'Bank Loan - FFR Sensitivity (se)&');
fprintf(FID, num2str(VAR{ismember(VAR(:,1),{'BankLoan-FFR 3-year sensitivity: '}),3},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, 'Repricing Duration&');
fprintf(FID, num2str(balance_sheet{ismember(balance_sheet(:,2),{'assets_dur_prepayment_adjusted'}),1},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, 'Repricing Duration (sd)&');
fprintf(FID, num2str(balance_sheet{ismember(balance_sheet(:,2),{'assets_dur_sd'}),1},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, 'Reserve Ratio&');
fprintf(FID, num2str(balance_sheet{ismember(balance_sheet(:,2),{'effective_reserve_ratio'}),1},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, 'Depositor Sensitivity to Deposit Rates&');
fprintf(FID, num2str(-demand_dep{ismember(demand_dep(:,1),{'alpha'}),2},'%.3f')); % flip the sign to change the sensitivity to spreasds to sensitivity to rate
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, 'Depositor Sensitivity to Deposit Rates (se)&');
fprintf(FID, num2str(demand_dep{ismember(demand_dep(:,1),{'alpha (se)'}),2},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, 'Dispersion of Depositor Sensitivity to Deposit Rates&');
fprintf(FID, num2str(demand_dep{ismember(demand_dep(:,1),{'sigma_alpha'}),2},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, 'Dispersion of Depositor Sensitivity to Deposit Rates (se)&');
fprintf(FID, num2str(demand_dep{ismember(demand_dep(:,1),{'sigma_alpha (se)'}),2},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, 'Borrower Sensitivity to Loan Rates&');
fprintf(FID, num2str(demand_loan{ismember(demand_loan(:,1),{'alpha'}),2},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, 'Borrower Sensitivity to Loan Rates (se)&');
fprintf(FID, num2str(demand_loan{ismember(demand_loan(:,1),{'alpha (se)'}),2},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
if i<=3
    fprintf(FID, 'Convenience of Holding Cash ($q_c$)&');
    fprintf(FID, num2str(demand_dep{ismember(demand_dep(:,1),{'q_cash'}),2},'%.3f'));
    fprintf(FID, '&');
    fprintf(FID, '\\\\\n');
    fprintf(FID, 'Convenience of Holding Cash (se)&');
    fprintf(FID, num2str(demand_dep{ismember(demand_dep(:,1),{'q_cash (se)'}),2},'%.3f'));
    fprintf(FID, '&');
    fprintf(FID, '\\\\\n');
    fprintf(FID, 'Convenience of Holding Deposits ($q_d$)&');
    fprintf(FID, num2str(demand_dep{ismember(demand_dep(:,1),{'q_deposit_sum'}),2},'%.3f'));
    fprintf(FID, '&');
    fprintf(FID, '\\\\\n');
    fprintf(FID, 'Convenience of Holding Deposits (se)&');
    fprintf(FID, num2str(demand_dep{ismember(demand_dep(:,1),{'q_deposit_sum (se)'}),2},'%.3f'));
%     fprintf(FID, '&');
%     fprintf(FID, '\\\\\n');
%     fprintf(FID, 'Convenience of Borrowing Bonds&');
%     fprintf(FID, num2str(demand_loan{ismember(demand_loan(:,1),{'q_NonBankLoan'}),2},'%.3f'));
%     fprintf(FID, '&');
%     fprintf(FID, '\\\\\n');
%     fprintf(FID, 'Convenience of Borrowing Bonds (se)&');
%     fprintf(FID, num2str(demand_loan{ismember(demand_loan(:,1),{'q_NonBankLoan (se)'}),2},'%.3f'));
    fprintf(FID, '&');
    fprintf(FID, '\\\\\n');
    fprintf(FID, 'Convenience of Borrowing Loans ($q_l$)&');
    fprintf(FID, num2str(demand_loan{ismember(demand_loan(:,1),{'q_BankLoan_sum'}),2},'%.3f'));
    fprintf(FID, '&');
    fprintf(FID, '\\\\\n');
    fprintf(FID, 'Convenience of Borrowing Loans (se)&');
    fprintf(FID, num2str(demand_loan{ismember(demand_loan(:,1),{'q_BankLoan_sum (se)'}),2},'%.3f'));
    fprintf(FID, '&');
    fprintf(FID, '\\\\\n');
else % use the whole sample outside option estimates for large and small banks
    fprintf(FID, 'Convenience of Holding Cash ($q_c$)&');
    fprintf(FID, num2str(quality_dep{ismember(quality_dep(:,1),{'q_cash'}),2},'%.3f'));
    fprintf(FID, '&');
    fprintf(FID, '\\\\\n');
    fprintf(FID, 'Convenience of Holding Cash (se)&');
    fprintf(FID, num2str(quality_dep{ismember(quality_dep(:,1),{'q_cash (se)'}),2},'%.3f'));
    fprintf(FID, '&');
    fprintf(FID, '\\\\\n');
    fprintf(FID, 'Convenience of Holding Deposits ($q_d$)&');
    fprintf(FID, num2str(quality_dep{ismember(quality_dep(:,1),{'q_deposit_sum'}),2},'%.3f'));
    fprintf(FID, '&');
    fprintf(FID, '\\\\\n');
    fprintf(FID, 'Convenience of Holding Deposits (se)&');
    fprintf(FID, num2str(quality_dep{ismember(quality_dep(:,1),{'q_deposit_sum (se)'}),2},'%.3f'));
    fprintf(FID, '&');
    fprintf(FID, '\\\\\n');
%     fprintf(FID, 'Convenience of Borrowing Bonds&');
%     fprintf(FID, num2str(quality_loan{ismember(quality_loan(:,1),{'q_NonBankLoan'}),2},'%.3f'));
%     fprintf(FID, '&');
%     fprintf(FID, '\\\\\n');
%     fprintf(FID, 'Convenience of Borrowing Bonds (se)&');
%     fprintf(FID, num2str(quality_loan{ismember(quality_loan(:,1),{'q_NonBankLoan (se)'}),2},'%.3f'));
%     fprintf(FID, '&');
%     fprintf(FID, '\\\\\n');
    fprintf(FID, 'Convenience of Borrowing Loans ($q_l$)&');
    fprintf(FID, num2str(quality_loan{ismember(quality_loan(:,1),{'q_BankLoan_sum'}),2},'%.3f'));
    fprintf(FID, '&');
    fprintf(FID, '\\\\\n');
    fprintf(FID, 'Convenience of Borrowing Loans (se)&');
    fprintf(FID, num2str(quality_loan{ismember(quality_loan(:,1),{'q_BankLoan_sum (se)'}),2},'%.3f'));
    fprintf(FID, '&');
    fprintf(FID, '\\\\\n');
end

fprintf(FID, 'Number of Banks&');
fprintf(FID, num2str(HHI{ismember(HHI(:,2),{'J_hat'}),1},'%.3f'));
fprintf(FID, '&');
fprintf(FID, '\\\\\n');
fprintf(FID, '\\hline ');
fprintf(FID, '\\hline ');
fprintf(FID, '\\end{tabular*}\n');
% fprintf(FID, strcat('\\caption{','Sample',num2str(i),'}'));
fclose(FID);

    
end
%fprintf(FID, '\\end{document}');







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demand Parameter Estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

load dep_1.mat
load loan_1.mat % load the demand estimates from loan market

Table_b=[];
Table_se=[];

Table_b=[ [-(beta_hat(1)+sum(qv.*qweight)*theta_hat(1)); beta_hat(2:nchar);theta_hat(1:end)/sqrt(12)]  ... 
    [btsls(1:nchar);zeros(nrc,1)] ];                                        %divide by sqrt(12) to convert the range of a uniform distribution to standard deviation                      
Table_se=[ [se12gmm(1:nchar);se12gmm(end-nrc+1:end)/sqrt(12)] ... 
    [setsls(1:1);zeros(nrc,1);setsls(2:nchar);]];       %divide by sqrt(12) to convert the range of a uniform distribution to standard deviation

[r c]=size(Table_se);
Table_demand=cell(3*r,c);

[r c]=size(Table_se);
Table_se2=cell(r,c);
Table_b2=cell(r,c);


for i=1:r
    for j=1:c
        zscore=abs(Table_b(i,j))/Table_se(i,j);
        if zscore<=1.644852;     star=''; end
        if zscore>1.644852;     star='*'; end
        if zscore>1.959964;     star='**'; end
        if zscore>2.575829;     star='***'; end
        Table_b2{i,j}=strcat(num2str(Table_b(i,j),'%0.3f'),star); %index of cell matrix needs {}
        Table_se2{i,j}=strcat('[',num2str(Table_se(i,j),'%0.3f'),']');
    end
end
Table_b2{4,2}='';
Table_se2{4,2}='';


var_name={'Yield Sensitivity ($\\alpha$)'  ...
    'Log No. of Branches ($\\beta_1$)' 'Log No. of Employees ($\\beta_2$)' ...
    'Yield Sensitivity Dispersion ($\\sigma_{\\alpha}$)' };

filepath = '.\intermediate_output\';
texfolder=strcat('../','output\Tables');
FID = fopen(strcat(texfolder,'\demand_blp.tex'), 'w');
fprintf(FID, '\\def\\sym#1{\\ifmmode^{#1}\\else\\(^{#1}\\)\\fi}');
fprintf(FID, '\\begin{tabular*}{\\hsize}{@{\\hskip\\tabcolsep\\extracolsep\\fill}l*{2}{c}}');
fprintf(FID, '\\hline ');
fprintf(FID, '\\hline ');
fprintf(FID, '&Deposit &Loan  \\\\ [1ex] \\hline  ');

for k=1:4 % this omit the time invariant characteristics
    fprintf(FID, var_name{k});
    fprintf(FID, '&');
    fprintf(FID, Table_b2{k,1});
    fprintf(FID, '&');
    fprintf(FID, Table_b2{k,2});
    fprintf(FID, '\\\\');
    fprintf(FID, '\n');

    
    fprintf(FID, '');
    fprintf(FID, '&');
    fprintf(FID, Table_se2{k,1});
    fprintf(FID, '&');
    fprintf(FID, Table_se2{k,2});
    fprintf(FID, '\\\\');
    fprintf(FID, '\n');

    fprintf(FID, '');
    fprintf(FID, '\\\\');
    fprintf(FID, '\n');
    
end
fprintf(FID, '\\hline ');
fprintf(FID, 'Sector F.E.&Y&Y \\\\ \n');
fprintf(FID, 'Time F.E.&Y&Y \\\\ \n');
fprintf(FID, 'Observations& %.0f &  %.0f',nobs',nobs');
fprintf(FID, '\\\\');
fprintf(FID, '\n');
fprintf(FID, 'Adj. Rsq& %.3f & %.3f',Rsq_gmm',Rsq_tsls');
fprintf(FID, '\\\\');
fprintf(FID, '\n');
fprintf(FID, '\\hline ');
fprintf(FID, '\\hline ');


fprintf(FID, '\\end{tabular*}\n');
fclose(FID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clean up folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% for i=1:5
%     filename=strcat(texfolder,'loan_',num2str(i),'.mat');
%     if exist(filename,'file')==2     
%         delete(filename);
%     end
% end
% 
% 
% for i=1:5
%     filename=strcat(texfolder,'deposit_',num2str(i),'.mat');
%     if exist(filename,'file')==2     
%         delete(filename);
%     end
% end
% 
% 
% filename='beta.mat';
% if exist(filename,'file')==2
%     delete(filename);
% end
% 
% 
% filename='mvalold.mat';
% if exist(filename,'file')==2
%     delete(filename);
% end