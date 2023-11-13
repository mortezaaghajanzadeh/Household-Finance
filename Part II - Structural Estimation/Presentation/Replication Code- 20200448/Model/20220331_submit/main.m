	%	First define parameters:		
diary off				
diary (strcat('output_',date,'.txt'))				
addpath('BLP');				
addpath('Library');				
addpath('Data');				
addpath('Results');	
addpath('Script');	

Par.cd0        = 0.0095;                                                       	%	The linear cost of keeping deposit
Par.cd1        = 0.000;                                                     	%	The quadratic cost of keeping deposit is Par.cd1/2*D^2
Par.cb0        = 0.0066;                                                       	%	The linear cost of generating loans
Par.cb1        = 0.000;                                                     	%	The quadratic cost of generating loans is is Par.cb1/2*B^2
				
Par.ND          = 6;                                                            %	Number of competing banks on deposit market
Par.NB          = 6;                                                            %	Number of competing banks on loan market

Par.W0       = 0.2167*Par.ND/Par.NB;                                            %   Initial wealth of households 
                                                                      				
		%	Estimates from BLP
				
Par.alphaD       = 96.8;                                                       	%	Mean Interest Rate elasticity on the deposit market
Par.sigalphaD  = 191.6;                                                         %	Dispersion of Interest Rate elasticity on the deposit market
Par.alpha         = 146.22;                                                   	%	Interest Rate elasticity on the loan market
				
Par.lD         = 1.985 + 1.455 - log(Par.ND);                                   %   Liquidity value of bank deposit on the deposit market
Par.lM        = 1.985;                                                      	%	Liquidity value of money on the deposit market
Par.lB         = 0.00;                                                          %	Liquidity value of bond
Par.os         = -9.631;                                                        %	Liquidity value of outside option
Par.nuB      = 1.1509 - log(Par.NB) - Par.os;                                   %	Net Return of bank debt on the loan market
Par.nuC      = 0 - Par.os;                                               	  	%	Net Return of coupon bond on the loan market
Par.e          = 0;			
Par.issue    = 0.01;							

		%	Shocks
Par.rhoA       = 0.6;                                                      	    %	Persistence of chargeoff
Par.nuA        = -1.3 - log(100);                                               %   Mean of banks' chargeoff
Par.sigA       = 1.2*sqrt(1-Par.rhoA^2);                                        %	Standard deviation of innnovation to the chargepff
			
Par.rhof       = 0.9000;                                                     	%	Persistence of fed funds rate
Par.nuf        = 0.3 - log(100);                                                %	Mean fed funds rate
Par.sigf       = sqrt(1.6*(1-Par.rhof^2));                                      %	Standard deviation of innnovation to the fed funds rate
Par.rhoAf0      = -0.11;                                                   		%	(\partial LOG chargeoff)/ (\partial LOG FFR)
Par.rhoAf1      = -0.04;                                                        %	(\partial chargeoff)/ (\partial FFR), ajudted for the scaling
		%	Bank's Problem
Par.beta       = 0.9545;                                                     	%	The bank's discount factor
Par.ksi          = 0.06;                                                     	%   Capital cosntraint
Par.theta      = 0.02495;                                                       %	Deposit constraint
Par.theta2    = 0.00;				
Par.mu         = 0.29;                                                     		 %	Fraction of maturing debt
Par.phi         = 1.00;                                                     	 %	Fraction of repricing loans
Par.a0          = 0.00;                                                     	 %	Linear cost of holding nonreservables in excess of the Fed funds rate
Par.a1          = 0.010;                                                     	 %	Quadratic cost of holding nonreservables (Par.a0*N+Par.a1/2*N^2)
Par.tauc       = 0.35;                                                     		 %	Corporate tax rate
Par.taup       = 0.20;				
Par.fix          = 0.8*1e-4;				
		
		%	Other parameters	
Par.nA       = 7;				
Par.nf        = 11;				
Par.nL        = 21;				
Par.nE       = 21;				
Par.nD       = 21;              				
Par.nL_Expand  = 21;              				
Par.nE_Expand  = 21;   				
			
Par.howard   = 0;                                                             %	   Counts in Howard Improment Method
Par.mcqueen= 0;                                                               %	   Whether to use Mcqueen Porteurs Bounds
Par.mA         = 1;                                                           %    How many additional grid points to be inserted in the simulation step
Par.mf          = 5;                                                          %    How many additional grid points to be inserted in the simulation step
Par.mL         = 1;				
Par.mE         = 1;				
			
Par = Grid_exp(Par);			
	
% plotyy(Par.fGrid,Par.D_star,Par.fGrid,Par.fGrid-Par.rd_star) 		



%% Now, solve the model and simulate on origional grid points

% Generate and store seeds for random shocks; this code only needs to be
% run once unless we change the dimensions for the simulation
% 
% N = 6;
% t = 50000;
% gen_seed(N,t);
tic
shock0= genshocks(Par); save('Data\shocks_6_500000_1_0_0','shock0');
solution0 = solve_model_MP(Par,'YES DMP','YES LMP', [], [], 'SAVE Kconst'); 

load('Data\shocks_6_500000_1_0_0'); 
[firm0, simu_con0, simu_un0, rec0] = SimulatePanel_short(Par, solution0, shock0);
elapsetime = toc;
str = ['Time to solve the model = ', num2str(elapsetime)]; disp(str)

%% Print the results

disp("Panel A: Calibrated Parameters")
disp([Par.tauc, Par.theta, Par.ksi, Par.NB])

disp("Panel B: Parameter Estimated Separately")
disp([1/Par.mu, Par.nuf, Par.sigf, Par.rhof, Par.nuA, Par.sigA, Par.rhoA, Par.rhoAf1*0.022/0.008]);

disp("Panel C: Parameter Estimated via BLP")
disp([Par.alphaD/100, Par.sigalphaD/100/sqrt(12), -Par.alpha/100, Par.lD+log(Par.ND), Par.lM, Par.nuB+Par.os+log(Par.NB)]);

disp("Panel D: Parameter Estimated via SMD")
disp([1/Par.beta-1, Par.W0, Par.os, Par.a1, Par.cd0, Par.cb0, Par.fix/simu_un0.meanE]); 

disp("Actual (first line) versus simulated (second line) moments")

mmt_actual = [0.0338, -0.299, 0.126, 0.0129, 0.0203, 0.699, 0.0120, 11.20, 2.061, -0.995, -1.592];
mmt_std = [  0.006  0.021  0.026  0.001  0.002  0.037  0.001  0.494  0.242  0.126  0.376];

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

disp([mmt_actual; mmt0; (-mmt_actual+mmt0)./mmt_std]')


%% Check R-sqr

mean0  = sum(Par.qv.*Par.qweight);
sensitivity0 = Par.alphaD + Par.sigalphaD*(Par.qv - mean0); sensitivity0 = sensitivity0 + (1-sensitivity0).*(sensitivity0<1);

imax = 10000;  % dont use the full simulation to save time
N = 6;
rd_bar1 = zeros(N*imax,Par.nq);
for i = 1:imax
    for j = 1:N
        for k = 1:Par.nq        
            indexk = find([1:N] ~= j);
            temp = exp(sensitivity0(k)*firm0.rd1(i,indexk) + Par.lD);
            exp_rate = sum(temp)/length(indexk); 
        	log_rate= (log(exp_rate) - Par.lD)/(sensitivity0(k));
            rd_bar1((i-1)*N+j,k)  = log_rate; 
        end
    end
end

rd_bar1 = reshape(rd_bar1,[],1);
rd_bar0 = solution0.rd_bar(:,shock0.indexf(1:imax,:)); 
rd_bar0 = kron(rd_bar0',ones(N,1)); rd_bar0 = reshape(rd_bar0,[],1);

stats = 1-(sqrt(sum((rd_bar1-rd_bar0).^2))/sqrt(sum((rd_bar1).^2)));
str = ['Deposit Mkt Forecast,  R-sqr = ', num2str(stats(1))]; disp(str)


r_bar1 = zeros(N*imax,1);
for i = 1:imax 
    for j = 1:N
        indexk = find([1:N] ~= j);
        temp = exp(Par.nuB - Par.alpha*reshape(firm0.r1(i,indexk),1,[]));
        exp_rate = sum(temp)/length(indexk);   
        log_rate = (log(exp_rate) - Par.nuB)/(- Par.alpha);
        r_bar1((i-1)*N+j) = log_rate;     
    end
end

r_bar1 = reshape(r_bar1,[],1);
r_bar0 = solution0.r_bar(:,shock0.indexf(1:imax,:)); 
r_bar0 = kron(r_bar0',ones(N,1)); r_bar0 = reshape(r_bar0,[],1);
stats = 1-(sqrt(sum((r_bar1-r_bar0).^2))/sqrt(sum((r_bar1).^2)));
str = ['Loan Mkt Forecast,  R-sqr = ', num2str(stats(1))]; disp(str)

save('Results\main')
