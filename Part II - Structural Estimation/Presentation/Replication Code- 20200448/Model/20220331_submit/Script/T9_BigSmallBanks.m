%% Run main.m and get the full sample results
clear all
% main
load('Results\main')  % can directly load solution; otherwise, run main.m


%% Big Banks

ParB = Par;
ParB.cd0        = 0.0103;                                                 % The linear cost of keeping deposit
ParB.cb0        = 0.0060;                                                 % The linear cost of generating loans

ParB.a0         = 0.00;                                                      % Linear cost of holding nonreservables in excess of the Fed funds rate
ParB.a1         = 0.006;                                                    % Quadratic cost of holding nonreservables (Par.a0*N+Par.a1/2*N^2)
ParB.fix        = 0.4*1e-4;

ParB.W0       = 0.228*Par.ND/Par.NB;                              %   Initial wealth of households 

ParB = Grid_exp(ParB);
solutionB = solve_model_MP(ParB,'YES DMP','YES LMP');
load('Data\shocks_6_500000_1_0_0'); [~, simu_conB, simu_unB] = SimulatePanel_short(ParB, solutionB, shock0);

% Small Banks

ParS = Par;
ParS.cd0        = 0.0077;                                                  % The linear cost of keeping deposit
ParS.cb0        = 0.0090;                                                  % The linear cost of generating loans

ParS.a0         = 0.000;                                                     % Linear cost of holding nonreservables in excess of the Fed funds rate
ParS.a1         = 0.015;                                                     % Quadratic cost of holding nonreservables (Par.a0*N+Par.a1/2*N^2)
ParS.fix        = 1.7*1e-4;

ParS.W0       = 0.140*Par.ND/Par.NB;                              %   Initial wealth of households 


ParS = Grid_exp(ParS);
solutionS = solve_model_MP(ParS,'YES DMP','YES LMP');
load('Data\shocks_6_500000_1_0_0'); [~, simu_conS, simu_unS] = SimulatePanel_short(ParS, solutionS, shock0);

% Intermediate case: close financing friction

ParI = ParB;
ParI.a1         = 0.015;                                                % Run the intermadiate case where all parameters are set equal to that among small banks, only change the financing costs

ParI = Grid_exp(ParI);
solutionI = solve_model_MP(ParI,'YES DMP','YES LMP');
load('Data\shocks_6_500000_1_0_0'); [~, simu_conI, simu_unI] = SimulatePanel_short(ParI, solutionI, shock0);

save('Results\T9_BigSmall','ParS', 'ParI', 'ParB','simu_unB', 'simu_conB','simu_unI', 'simu_conI','simu_unS', 'simu_conS')


%% Report the moments

load('Results\T9_BigSmall')

mmtBig = [0.0360, -0.355, 0.133, 0.0132, 0.0177, 0.666, 0.0096, 11.36, 2.06, -0.995];
mmtSmall =  [0, -0.153, 0.087, 0.0125, 0.02716, 0.784, 0.0190, 10.78, 0, -0.995];
mmtBig_std = [0.002  0.019  0.025  0.001  0.001  0.032  0.001  0.580  0.077  0.299];
mmtSmall_std = [0.002  0.005  0.009  0.001  0.001  0.012  0.001  0.198  0.086  0.299];

mmtB = [simu_unB.meanC2V, simu_unB.meanN2D, simu_unB.stdN2D, simu_unB.spreadrd, simu_unB.spreadr, simu_unB.meanD2A, ...
              simu_unB.meanFIX1, simu_unB.meanLEV, simu_unB.meanM2B, simu_unB.reg_T_agg];
mmtS = [simu_unS.meanC2V*0, simu_unS.meanN2D, simu_unS.stdN2D, simu_unS.spreadrd, simu_unS.spreadr, simu_unS.meanD2A, ...
              simu_unS.meanFIX1, simu_unS.meanLEV, simu_unS.meanM2B*0, simu_unS.reg_T_agg];   % we dont have "dividend payout" and "market to book ratio" for small banks beucase they are mostly private
          
disp('Moment condition for big and small banks:')
disp([mmtBig', mmtB', (-mmtBig'+mmtB')./mmtBig_std', mmtSmall', mmtS', (-mmtSmall'+mmtS')./mmtSmall_std'])


%% Report the parameters 

str = (['Big banks par = ', num2str([ParB.a1, ParB.cd0, ParB.cb0, ParB.fix/simu_unB.meanE, ParB.W0]), ', Sensitivity of loans = ', num2str(simu_unB.reg_B_agg)]); disp(str)
str = (['Small banks par = ', num2str([ParS.a1, ParS.cd0, ParS.cb0, ParS.fix/simu_unS.meanE, ParS.W0]), ', Sensitivity of loans = ', num2str(simu_unS.reg_B_agg)]); disp(str)

percentage_financing = (simu_unI.reg_B_agg - simu_unB.reg_B_agg)/(simu_unS.reg_B_agg - simu_unB.reg_B_agg);
str = (['If we close external financing gap, the loan sensitivity will become ', num2str(simu_unI.reg_B_agg)]); disp(str)
str = (['The external financing cost accounts for ', num2str(percentage_financing), ' of the loan sensitivity dispersion']); disp(str)
