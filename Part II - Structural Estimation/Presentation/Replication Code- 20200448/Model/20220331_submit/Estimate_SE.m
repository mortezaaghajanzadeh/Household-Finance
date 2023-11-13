%% Select the list of parameters and moments

clear all
main	

for i  = 1:5
subsample = i;
flag = 0;

if subsample ==1 % Full Sample
    load('Results\main')  % can directly load solution; otherwise, run main.m  
    mmt_actual = [0.0338, -0.299, 0.126, 0.0129, 0.0203, 0.699, 0.0120, 11.20, 2.061, -0.995, -1.592];
    addpath('Data');
    addpath('Library');
    [ID, infl, infl_count] = getinfl('standard_error_1.csv');
    
elseif subsample ==2 % Early
    load('Results\T10_EarlyLate')
    mmt_actual = [0.0310, -0.342, 0.103, 0.0195, 0.0183, 0.679, 0.0130, 12.46, 2.767, -0.995, 0];
    Par = Par_Early;
    simu_un0 = simu_unE;
    addpath('Data');
    addpath('Library');
    [ID, infl, infl_count] = getinfl('standard_error_2.csv');
    
elseif subsample ==3 % Late
    load('Results\T10_EarlyLate')
    mmt_actual = [0.0360, -0.256, 0.106, 0.0064, 0.0224, 0.719, 0.0120, 9.933, 1.538, -0.995, 0];
    Par = Par_Late;
    simu_un0 = simu_unL;
    addpath('Data');
    addpath('Library');
    [ID, infl, infl_count] = getinfl('standard_error_3.csv');
    
elseif subsample ==4 % Big
    load('Results\T9_BigSmall')
    mmt_actual = [0.0360, -0.355, 0.133, 0.0132, 0.0267, 0.666, 0.0096, 11.36, 2.06, -0.995, 0];
    Par = ParB;
    simu_un0 = simu_unB;
    addpath('Data');
    addpath('Library');
    [ID, infl, infl_count] = getinfl('standard_error_4.csv');
    
elseif subsample ==5 % Small
    load('Results\T9_BigSmall')
    mmt_actual = [0, -0.153, 0.087, 0.0125, 0.0319, 0.784, 0.0190, 10.78, 0, -0.995, 0];
    Par = ParS;
    simu_un0 = simu_unS;
    addpath('Data');
    addpath('Library');
    [ID, infl, infl_count] = getinfl('standard_error_5.csv');
end

mmt0= [simu_un0.meanC2V, ...
             simu_un0.meanN2D, ...
             simu_un0.stdN2D, ...
             simu_un0.spreadrd, ...
             simu_un0.spreadr, ...
             simu_un0.meanD2A, ...
             simu_un0.meanFIX1, ...
             simu_un0.meanLEV, ...
             simu_un0.meanM2B, ...
             simu_un0.reg_T_agg, ...
             simu_un0.reg_B_agg];

Nmmt = length(mmt0);    

if subsample <= 3
    i_grid = [1,2,3,4,5,6,7]; % We estimate 7 parameters in Early v.s. Late subsamples
    Npar = length(i_grid);
else 
    i_grid = [2,4,5,6,7];          % We only re-estimate the bank operatoinal and financing cost parameters in Bag v.s. Small subsamples
    Npar = length(i_grid);
end
    
% 1. Par.beta
% 2. Par.W0
% 3. Par.os
% 4. Par.a1
% 5. Par.cb0
% 6. Par.cd0
% 7. Par.fix


% Construct influence fucntion

count_matrix = repmat(infl_count, 9,1).*repmat(infl_count',1,9);
plain_covM   = cov(infl)./sqrt(count_matrix);
plain_covM(10,10) = 0.299^2;     % the total credit-FFR sensitivity and the bank loan-FFR sensitivity are calculated based on aggregate (not firm level) data
plain_covM(11,11) = 0.503^2;     % so we directly take the s.e. estimates from the VAR (and not use influence functions on those).


% get demean and clustered weight matrix

[demean_infl, cluster_infl] = demean(ID, infl);
Nob = size(infl,1);
Nob_cluster = size(cluster_infl,1);

cluster_covM   = cov(cluster_infl)./sqrt(count_matrix*Nob_cluster/Nob);
demean_covM = cov(demean_infl)./sqrt(count_matrix);
cluster_covM(10,10)= 0.299^2;
demean_covM(10,10) = 0.299^2;

cluster_covM(11,11)= 0.503^2;
demean_covM(11,11) = 0.503^2;

W0 = inv(cluster_covM);



% Estimate parameters

% if set flag == 1 then we estimate the parameters 
% if set flag == 0 (default) then we take the supplied parameters and calculate the standard errors 

s0 = [Par.beta, Par.W0, Par.os, Par.a1, Par.cb0, Par.cd0, Par.fix];
s_disp = [1/Par.beta-1, Par.W0, Par.os, Par.a1, Par.cb0, Par.cd0, Par.fix/simu_un0.meanE];

x0  = s0; 
myfun=@(x) calculate_err(Par,x,W0,mmt_actual); 

if flag == 1
    options = optimset('Display','iter','PlotFcns',@optimplotfval);
    x0  = fminsearch(myfun,x0,options);
end
[~,mmtx]=myfun(x0);


%  Estimate the VarCov matrix for parameters and moments

delta_grid =  [0.01 0.025 0.05 0.10];  % different "steps" while taking numrical derivatives
se_par_matrix = zeros(length(delta_grid),Npar);
se_mmt_matrix= zeros(length(delta_grid),Nmmt);

for k = 1:length(delta_grid)
   
    delta = delta_grid(k);
    jacob = zeros(Nmmt,Npar);
    
    for iindex = 1:Npar

        i = i_grid(iindex);
        Pari = Par;
        if i ==1
            Pari.beta = Par.beta*(1-delta);  deltai = delta*Par.beta;
        elseif i ==2
            Pari.W0 = Par.W0*(1-delta);     deltai = delta*Par.W0;
        elseif i ==3
            Pari.os = Par.os*(1-delta*10);  deltai = delta*10*Par.os;
            Pari.nuB  = 1.1509 - log(Par.NB) - Pari.os;                                	 
            Pari.nuC  = 0 - Pari.os;                                               	
        elseif i ==4
            Pari.a1 = Par.a1*(1-delta);        deltai = delta*Par.a1;
        elseif i ==5
            Pari.cb0 = Par.cb0*(1-delta);    deltai = delta*Par.cb0;          
        elseif i ==6
            Pari.cd0 = Par.cd0*(1-delta);    deltai = delta*Par.cd0; 
        elseif i ==7
            Pari.fix = Par.fix*(1-delta);        deltai = delta*Par.fix; 
        end       
        
        solutioni = solve_model_MP(Pari,'YES DMP','YES LMP', [], [], 'SAVE Kconst');
        [firmi, simu_coni, simu_uni, reci] = SimulatePanel_short(Pari, solutioni, shock0);
       
        mmti= [simu_uni.meanC2V, ...
             simu_uni.meanN2D, ...
             simu_uni.stdN2D, ...
             simu_uni.spreadrd, ...
             simu_uni.spreadr, ...
             simu_uni.meanD2A, ...
             simu_uni.meanFIX1, ...
             simu_uni.meanLEV, ...
             simu_uni.meanM2B, ...
             simu_uni.reg_T_agg, ...
             simu_uni.reg_B_agg];

        jacob(:,iindex) = reshape(mmt0-mmti,[],1)/deltai;
    end     
    sandwich_side = (jacob'*W0*jacob)\jacob'*W0;
    
    cov_par   = sandwich_side*cluster_covM*sandwich_side'; % small "sandwich" to calculate the s.e. for parameters
    cov_mmt = (eye(Nmmt)-jacob*sandwich_side)*cluster_covM*(eye(Nmmt)-jacob*sandwich_side)'; % big "sandwich" to calculate the s.e. for moments 
                                                                                                                                                                     % following Liu, Whited, Zhang (2009, JPE)
    se_par_matrix(k,:) = sqrt(reshape(diag(cov_par),1,[]));
    se_mmt_matrix(k,:) = sqrt(reshape(diag(cov_mmt),1,[]));
end

se_simple = sqrt(sum(eye(length(cluster_covM)).*(cluster_covM)));
se_mmt_matrix = [se_simple; se_mmt_matrix];

% save results
if subsample == 1
    save('Results\se_full','se_par_matrix','se_mmt_matrix','cov_par','cov_mmt','*_covM','s0','x0','mmt0','mmtx','mmt_actual')
elseif subsample ==2
    save('Results\se_early','se_par_matrix','se_mmt_matrix','cov_par','cov_mmt','*_covM','s0','x0','mmt0','mmtx','mmt_actual')
elseif subsample ==3
    save('Results\se_late','se_par_matrix','se_mmt_matrix','cov_par','cov_mmt','*_covM','s0','x0','mmt0','mmtx','mmt_actual')
elseif subsample ==4
    save('Results\se_big','se_par_matrix','se_mmt_matrix','cov_par','cov_mmt','*_covM','s0','x0','mmt0','mmtx','mmt_actual')
elseif subsample ==5
    save('Results\se_small','se_par_matrix','se_mmt_matrix','cov_par','cov_mmt','*_covM','s0','x0','mmt0','mmtx','mmt_actual')
end

end

%% Generate outputs for the full sample
for indexi = 1:4

clearvars -except indexi

load('Results\se_full')
load('Results\main')
se_par = se_par_matrix(indexi,[1,2,3,4,6,5,7]);  % change the sequence to be consistent with the ordering oreporetd in the table 
se_par(7) = se_par(7)/simu_un0.meanE;          % scale
    
fileID = fopen('Results\standard error.txt', 'a');
fprintf(fileID,'\n\n Full sample parameters \n');
fprintf(fileID,'%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n ', se_par);
fclose(fileID);

se_mmt = se_mmt_matrix(5,:); 
fileID = fopen('Results\standard error.txt', 'a');
fprintf(fileID,'\n\n Full sample moment conditions \n');
fprintf(fileID,'%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n ', se_mmt);
fclose(fileID);

% Generate outputs for big and small banks

clearvars -except indexi

load('Results\se_big')
load('Results\T9_BigSmall')
se_par = se_par_matrix(indexi,[2,4,3,5,1]);  % change the sequence to be consistent with the ordering oreporetd in the table
se_par(4) = se_par(4)/simu_unB.meanE;          % scale
    
fileID = fopen('Results\standard error.txt', 'a');
fprintf(fileID,'\n\n Big bank parameters \n');
fprintf(fileID,'%6.3f %6.3f %6.3f %6.3f  %6.3f \n ', se_par);
fclose(fileID);

se_mmt = se_mmt_matrix(5,:); 
fileID = fopen('Results\standard error.txt', 'a');
fprintf(fileID,'\n\n Big bank moment conditions \n');
fprintf(fileID,'%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n ', se_mmt(1:end-1)); % We do not match the bank-lending FFR sensitivity in this case
fclose(fileID);                                                                                                                                          % instead, we use model to make predictions


load('Results\se_small')
load('Results\T9_BigSmall')
se_par = se_par_matrix(indexi,[2,4,3,5,1]);  % change the sequence to be consistent with the ordering oreporetd in the table
se_par(4) = se_par(4)/simu_unS.meanE;          % scale
    
fileID = fopen('Results\standard error.txt', 'a');
fprintf(fileID,'\n\n Small bank parameters \n');
fprintf(fileID,'%6.3f %6.3f %6.3f %6.3f  %6.3f \n ', se_par);
fclose(fileID);

se_mmt = se_mmt_matrix(5,:); 
fileID = fopen('Results\standard error.txt', 'a');
fprintf(fileID,'\n\n Small bank moment conditions \n');
fprintf(fileID,'%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n ', se_mmt(1:end-1));
fclose(fileID);


% Generate outputs for early and late subsamples

clearvars -except indexi

load('Results\se_early')
load('Results\T10_EarlyLate')
se_par = se_par_matrix(indexi,[1,2,4,6,5,7]);  % change the sequence to be consistent with the ordering oreporetd in the table 
se_par(6) = se_par(6)/simu_unE.meanE;          % scale
    
fileID = fopen('Results\standard error.txt', 'a');
fprintf(fileID,'\n\n Early subsample parameters \n');
fprintf(fileID,'%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n ', se_par);
fclose(fileID);

se_mmt = se_mmt_matrix(5,:); 
fileID = fopen('Results\standard error.txt', 'a');
fprintf(fileID,'\n\n  Early subsample moment conditions \n');
fprintf(fileID,'%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n ', se_mmt);
fclose(fileID);
                

clearvars -except indexi

load('Results\se_late')
load('Results\T10_EarlyLate')
se_par = se_par_matrix(indexi,[1,2,4,6,5,7]);  % change the sequence to be consistent with the ordering oreporetd in the table 
se_par(6) = se_par(6)/simu_unL.meanE;          % scale


fileID = fopen('Results\standard error.txt', 'a');
fprintf(fileID,'\n\n Late subsample parameters \n');
fprintf(fileID,'%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n ', se_par);
fclose(fileID);

se_mmt = se_mmt_matrix(5,:); 
fileID = fopen('Results\standard error.txt', 'a');
fprintf(fileID,'\n\n  Late subsample moment conditions \n');
fprintf(fileID,'%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n ', se_mmt);
fclose(fileID);
end
                
%     