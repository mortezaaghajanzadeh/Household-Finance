addpath('BLP');				
addpath('Library');				
addpath('Data');				
addpath('Results');	
addpath('Script');		


% First define parameters:		
diary off		
diary output_09oct2017.txt		
addpath('BLP');		
addpath('Library');		
addpath('Data');		
		
Par.fix          = 2.4*1e-4;				
Par.cd0        = 0.85/100;                                           
Par.cd1        = 0.000;                                                
Par.cb0        = 0.45/100;                                          
Par.cb1        = 0.000;                                                
Par.a0          = 0.000;                                               
Par.a1          = 0.010;                                               


%	
Par.ND               = 7;
Par.NB               = 7;                                                                                                          	
Par.alphaD      = 74.3;                                              
Par.sigalphaD = 146.7;                                             
Par.alpha        = 101.7;                                            

Par.lD             = 3.465 - log(Par.ND);                            
Par.lM            = 2.763;                                              
Par.lB              = 0.00;                                              
Par.os             = -9.631;                                                
Par.nuB          = -0.016 - log(Par.NB) - Par.os;          
Par.nuC          = 0 - Par.os;                                        
Par.e               = 0;		
Par.issue        = 0.01;		


Par.W0         = 0.1843;                                                  	
Par.mu         = 1/3.17;                                                     	

Par.beta       = 0.955;                                                     	
Par.ksi          = 0.06;                                                     	
Par.theta     = 0.028;                                                   	
Par.theta2   = 0.00;		
Par.phi         = 1.00;                                                     	
Par.tauc       = 0.35;                                                     
Par.taup      = 0.20;		


Par.rhoA       = 0.6;                                                      	
Par.nuA        = -1.4 - log(100);                                          	
Par.sigA        = 1.0*sqrt(1-Par.rhoA^2);                                                      	

Par.rhof       = 0.70;                                                     	
Par.nuf        = 1.3 - log(100);                                          	
Par.sigf        = sqrt(0.6*(1-Par.rhof^2));                                 	
Par.rhoAf0  = 0.04;                                                   	
Par.rhoAf1  = 0.01;                                                     	

	%	 Other parameters
	%	
Par.nA         = 7;		
Par.nf         = 11;		
Par.nL         = 21;		
Par.nE         = 21;		
Par.nD         = 21;              		
Par.nL_Expand  = 21;              		
Par.nE_Expand  = 21;   		
		
Par.howard     = 0;                                                       
Par.mcqueen  = 0;                                                         	
Par.mA        = 1;                                                         	
Par.mf         = 5;		
Par.mL         = 1;		
Par.mE         = 1;		
		
Par = Grid_adjust(Par);		
Par_Early = Par;

shock0= genshocks(Par_Early); save('Data\shocks_6_500000_1_0_0','shock0');
solutionE = solve_model_MP(Par_Early,'YES DMP','YES LMP', [], [], 'SAVE Kconst'); 
[firmE, simu_conE, simu_unE, recE] = SimulatePanel_short(Par_Early, solutionE, shock0);
mmt_E_actual = [0.0310, -0.342, 0.103, 0.0195, 0.0263, 0.679, 0.0130, 12.46, 2.767, -0.995, 0];

%% First, solve the model using parameters for the "Late" period

Par.fix          = 0.4*1e-4;				
Par.cd0        = 0.9/100;                                           
Par.cd1        = 0.000;                                                
Par.cb0        = 0.80/100;                                          
Par.cb1        = 0.000;                                                
Par.a0          = 0.000;                                               
Par.a1          = 0.010;                                                
                                           
Par.ND               = 5; 
Par.NB               = 5;                                                                                                          	

Par.alphaD      = 92.5;                                              
Par.sigalphaD  = 183;                                             
Par.alpha        = 145.4;                                            

Par.lD             = 2.34 - log(Par.ND);                            
Par.lM            = -0.444;                                              
Par.lB              = 0.00;                                              
Par.os             = -9.631;                                                 
Par.nuB          = 1.804 - log(Par.NB) - Par.os;          
Par.nuC          = 0 - Par.os;                                        
Par.e               = 0;		
Par.issue        = 0.01;			

Par.W0         = 0.2537;                                                  	
Par.mu         = 1/3.590;                                                     	

Par.beta       = 0.958;                                                     	
Par.ksi          = 0.06;                                                     	
Par.theta     = 0.022;                                                   	
Par.theta2   = 0.00;		
Par.phi         = 1.00;                                                     	
Par.tauc       = 0.35;                                                     
Par.taup      = 0.20;		

Par.rhoA       = 0.6;                                                      	
Par.nuA        = -1.2 - log(100);                                          	
Par.sigA       = 1.3*sqrt(1-Par.rhoA^2);                                                      	
		
Par.rhof       = 0.9;                                                     
Par.nuf        = -1.0- log(100);                                          	
Par.sigf       = sqrt(1.4*(1-Par.rhof^2));                                 	
Par.rhoAf0      = -0.16;                                                  	
Par.rhoAf1      = -0.20;                   

Par.nA         = 7;		
Par.nf         = 11;		
Par.nL         = 21;		
Par.nE         = 21;		
Par.nD         = 21;              		
Par.nL_Expand  = 21;              		
Par.nE_Expand  = 21;   		
		
Par.howard     = 0;                                                       
Par.mcqueen  = 0;                                                         	
Par.mA        = 1;                                                         	
Par.mf         = 5;		
Par.mL         = 1;		
Par.mE         = 1;		
		
Par = Grid_adjust(Par);		
Par_Late = Par;

shock0= genshocks(Par_Late);
solutionL = solve_model_MP(Par,'YES DMP','YES LMP', [], [], 'SAVE Kconst');
[firmL, simu_conL, simu_unL, recL] = SimulatePanel_short(Par_Late, solutionL, shock0);
mmt_L_actual = [0.0360, -0.256, 0.106, 0.0064, 0.0304, 0.719, 0.0120, 9.933, 1.538, -0.995, 0];



%% Group paramters into three categories and examine their effect on monetary transmission 

for i  = [1,2,3]
    
    if i == 1
        Par = Par_Late;
        Par.fix        = Par_Early.fix;				
        Par.cd0        = Par_Early.cd0;                                         
        Par.cb0        = Par_Early.cb0;                                  
        Par.a1         = Par_Early.a1;    
        Par.beta       = Par_Early.beta;                                                    	

    elseif i == 2  
        Par = Par_Late;
        Par.ND          = Par_Early.ND;        
        Par.NB          = Par_Early.NB;
        
        Par.alphaD      = Par_Early.alphaD;                                              
        Par.sigalphaD   = Par_Early.sigalphaD;                                       
        Par.alpha       = Par_Early.alpha;           
        
    elseif i ==3
        Par = Par_Late;
        Par.rhof       = Par_Early.rhof;                                                 	
        Par.nuf        = Par_Early.nuf;                                          	
        Par.sigf       = Par_Early.sigf;      
        
        Par.rhoAf0     = Par_Early.rhoAf0;                                                 	
        Par.rhoAf1     = Par_Early.rhoAf1; 
        
        Par.W0         = Par_Early.W0;                                            	
        Par.mu         = Par_Early.mu;                                                   	
        Par.rhoA       = Par_Early.rhoA;                                                     	
        Par.nuA        = Par_Early.nuA;                                       	
        Par.sigA       = Par_Early.sigA;    
        
        Par.lD         = Par_Early.lD;                            
        Par.lM         = Par_Early.lM;                                              
        Par.nuB        = Par_Early.nuB;          
        Par.nuC        = Par_Early.nuC;      
        
        Par.ksi        = Par_Early.ksi;                                                     	
        Par.theta      = Par_Early.theta;
    end                                       
                                                                                       	  
    Par = Grid_adjust(Par);		
        
    shock0= genshocks(Par); save('Data\shocks_6_500000_1_0_0','shock0');
    solution0 = solve_model_MP(Par,'YES DMP','YES LMP', [], [], 'SAVE Kconst');
    [firm0, simu_con0, simu_un0, rec0] = SimulatePanel_short(Par, solution0, shock0);
   
    loan_sensitivity(i) = rec0(4);
end

save('Results\T10_EarlyLate','Par_Early', 'Par_Late', 'simu_unE', 'simu_unL', 'loan_sensitivity')


%%   Print Results

load('Results\T10_EarlyLate')

disp("Panel A: Calibrated Parameters")
disp(['Early = ',num2str([Par_Early.tauc, Par_Early.theta, Par_Early.ksi, Par_Early.NB])])
disp(['Late = ',num2str([Par_Late.tauc, Par_Late.theta, Par_Late.ksi, Par_Late.NB])])
disp([' '])
disp("Panel B: Parameter Estimated Separately")
disp(['Early = ',num2str([1/Par_Early.mu, Par_Early.nuf, Par_Early.sigf, Par_Early.rhof, Par_Early.nuA, Par_Early.sigA, Par_Early.rhoA, Par_Early.rhoAf0])])
disp(['Late = ',num2str([1/Par_Late.mu, Par_Late.nuf, Par_Late.sigf, Par_Late.rhof, Par_Late.nuA, Par_Late.sigA, Par_Late.rhoA, Par_Late.rhoAf0])])
disp([' '])
disp("Panel C: Parameter Estimated via BLP")
disp(['Early = ',num2str([Par_Early.alphaD/100, Par_Early.sigalphaD/100/sqrt(12), -Par_Early.alpha/100, Par_Early.lD+log(Par_Early.ND), Par_Early.lM, Par_Early.nuB+Par_Early.os+log(Par_Early.NB)])])
disp(['Late = ',num2str([Par_Late.alphaD/100, Par_Late.sigalphaD/100/sqrt(12), -Par_Late.alpha/100, Par_Late.lD+log(Par_Late.ND), Par_Late.lM, Par_Late.nuB+Par_Late.os+log(Par_Late.NB)])])
disp([' '])
disp("Panel D: Parameter Estimated via SMD")
disp(['Early = ',num2str([1/Par_Early.beta-1, Par_Early.W0, Par_Early.os, Par_Early.a1, Par_Early.cd0, Par_Early.cb0, Par_Early.fix/simu_unE.meanE])])
disp(['Late = ',num2str([1/Par_Late.beta-1, Par_Late.W0, Par_Late.os, Par_Late.a1, Par_Late.cd0, Par_Late.cb0, Par_Late.fix/simu_unL.meanE])])
disp([' '])
                        
P1 = (loan_sensitivity(1)-simu_unL.reg_B_agg)/((loan_sensitivity(1)-simu_unL.reg_B_agg)+(loan_sensitivity(2)-simu_unL.reg_B_agg)+(loan_sensitivity(3)-simu_unL.reg_B_agg));
P2 = (loan_sensitivity(2)-simu_unL.reg_B_agg)/((loan_sensitivity(1)-simu_unL.reg_B_agg)+(loan_sensitivity(2)-simu_unL.reg_B_agg)+(loan_sensitivity(3)-simu_unL.reg_B_agg));
P3 = (loan_sensitivity(3)-simu_unL.reg_B_agg)/((loan_sensitivity(1)-simu_unL.reg_B_agg)+(loan_sensitivity(2)-simu_unL.reg_B_agg)+(loan_sensitivity(3)-simu_unL.reg_B_agg));

str = ['Contribution of operating costs = ', num2str(P1)]; disp(str)
str = ['Contribution of market power = ', num2str(P2)]; disp(str)
str = ['Contribution of low interest rate = ', num2str(P3)]; disp(str)