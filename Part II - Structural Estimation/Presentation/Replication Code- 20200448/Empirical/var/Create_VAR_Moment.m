% This example shows how to compute IRFs, HDs, and FEVDs in a VAR with 
% data for inflation, unemployment, and interest rates.  Identification 
% is achieved by imposing short-run restrictions, computed with a Cholesky 
% decomposition of the reduced-form residuals' covariance matrix.  

% The VAR Toolbox 2.0 is required to run this code. To get the 
% latest version of the toolboxes visit: 
% 
%       https://sites.google.com/site/ambropo/MatlabCodes
% 
% =======================================================================
% Ambrogio Cesa Bianchi, March 2015
% ambrogio.cesabianchi@gmail.com

clear all; clear session; close all; clc
warning off all

%% 1. PRELIMINARIES
% =======================================================================
% Add that folder plus all subfolders to the path.
addpath(genpath('.\'));
rng(0) 
% Load
[xlsdata, xlstext] = xlsread('Bernanke_Bliner_1959_2017.xlsx','Sheet1');

% Define data
X = xlsdata;
% Define label for plots
dates = xlstext(2:end,1);
vnames = xlstext(1,2:end);
% Define number of variables and of observations
[nobs, nvar] = size(X);


%% VAR ESTIMATION
% =======================================================================
% Set the case for the VARout (0, 1, or 2)
det = 0;
% Set number of nlags
nlags = 6;
% Estimate 
[VAR, VARopt] = VARmodel(X,nlags,det);
% Print at screen and create table
VARopt.vnames = vnames;
[beta, tstat, TABLE] = VARprint(VAR,VARopt);


%% IMPULSE RESPONSE
% =======================================================================
% Set options some options for IRF calculation
VARopt.nsteps = 48;
VARopt.ident = 'oir';
VARopt.quality = 0;
VARopt.pick = nvar;
% Compute IRF
[IRF, VAR] = VARir(VAR,VARopt);
% Compute error bands
[IRFINF,IRFSUP,IRFMED] = VARirband(VAR,VARopt);
% Plot
VARirplot(IRFMED,VARopt,IRFINF,IRFSUP);

disp(['TotalCredit-FFR 1-year sensitivity: ' num2str(IRFMED(12,1,nvar)/VAR.sigma(nvar,nvar)^.5*100) '%'])
disp(['BankLoan-FFR 1-year sensitivity: ' num2str(IRFMED(12,2,nvar)/VAR.sigma(nvar,nvar)^.5*100) '%'])
disp(['TotalCredit-FFR 2-year sensitivity: ' num2str(IRFMED(24,1,nvar)/VAR.sigma(nvar,nvar)^.5*100) '%'])
disp(['BankLoan-FFR 2-year sensitivity: ' num2str(IRFMED(24,2,nvar)/VAR.sigma(nvar,nvar)^.5*100) '%'])
disp(['TotalCredit-FFR 3-year sensitivity: ' num2str(IRFMED(36,1,nvar)/VAR.sigma(nvar,nvar)^.5*100) '%'])
disp(['BankLoan-FFR 3-year sensitivity: ' num2str(IRFMED(36,2,nvar)/VAR.sigma(nvar,nvar)^.5*100) '%'])

parameters={};
parameters{1,1}='TotalCredit-FFR 1-year sensitivity: ';parameters{1,2}=num2str(IRFMED(12,1,nvar)/VAR.sigma(nvar,nvar)^.5*100);
parameters{1,3}=num2str((IRFSUP(12,1,nvar)-IRFMED(12,1,nvar))/VAR.sigma(nvar,nvar)^.5*100/1.96); %se=(95_ci_up-mean)/1.96
parameters{2,1}='BankLoan-FFR 1-year sensitivity:';parameters{2,2}=num2str(IRFMED(12,2,nvar)/VAR.sigma(nvar,nvar)^.5*100);
parameters{2,3}=num2str((IRFSUP(12,2,nvar)-IRFMED(12,2,nvar))/VAR.sigma(nvar,nvar)^.5*100/1.96); %se=(95_ci_up-mean)/1.96
parameters{3,1}='TotalCredit-FFR 2-year sensitivity: ';parameters{3,2}=num2str(IRFMED(24,1,nvar)/VAR.sigma(nvar,nvar)^.5*100);
parameters{3,3}=num2str((IRFSUP(24,1,nvar)-IRFMED(24,1,nvar))/VAR.sigma(nvar,nvar)^.5*100/1.96); %se=(95_ci_up-mean)/1.96
parameters{4,1}='BankLoan-FFR 2-year sensitivity: ';parameters{4,2}=num2str(IRFMED(24,2,nvar)/VAR.sigma(nvar,nvar)^.5*100);
parameters{4,3}=num2str((IRFSUP(24,2,nvar)-IRFMED(24,2,nvar))/VAR.sigma(nvar,nvar)^.5*100/1.96); %se=(95_ci_up-mean)/1.96
parameters{5,1}='TotalCredit-FFR 3-year sensitivity: ';parameters{5,2}=num2str(IRFMED(36,1,nvar)/VAR.sigma(nvar,nvar)^.5*100);
parameters{5,3}=num2str((IRFSUP(36,1,nvar)-IRFMED(36,1,nvar))/VAR.sigma(nvar,nvar)^.5*100/1.96); %se=(95_ci_up-mean)/1.96
parameters{6,1}='BankLoan-FFR 3-year sensitivity: ';parameters{6,2}=num2str(IRFMED(36,2,nvar)/VAR.sigma(nvar,nvar)^.5*100);
parameters{6,3}=num2str((IRFSUP(36,2,nvar)-IRFMED(36,2,nvar))/VAR.sigma(nvar,nvar)^.5*100/1.96); %se=(95_ci_up-mean)/1.96

output_table_name=strcat('..\blp\intermediate_output\VAR','.xlsx');
xlswrite(output_table_name, parameters)

