function [Par_expand,solution_expand] = expand_solution(Par,solution0)

    % This function takes the original model solution defined on Par. It uses linear
    % interpolation to define it on a finer set of grids 
    % The finer set of grids are constructed by inserting (mA-1), (mf-1), (mE-1), (mL-1)
    % points in between each 2 points in the original grid space
    % The new model solution defined on the finer grid space is saved in the strut variable "solution_expand"


   
    % Define new space parameters #########################################
    
    Par_expand = Par;

    Par_expand.nA         = (Par.nA-1)*Par.mA+1;
    Par_expand.nf         = (Par.nf-1)*Par.mf+1;
    Par_expand.nL         = (Par.nL-1)*Par.mL+1;
    Par_expand.nE         = (Par.nE-1)*Par.mE+1;
    Par_expand.nL_Expand  = (Par.nL_Expand-1)*Par.mL+1;              
    Par_expand.nE_Expand  = (Par.nE_Expand-1)*Par.mE+1;   

    Par_expand = Grid_exp(Par_expand);


  
    % Interpolate #########################################################
 
    solution_expand.A0 = LinearInterp(Par, Par_expand, solution0.A0); 
    solution_expand.f0 = LinearInterp(Par, Par_expand, solution0.f0); 
    solution_expand.Ef0= LinearInterp(Par, Par_expand, solution0.Ef0); 
    solution_expand.L0 = LinearInterp(Par, Par_expand, solution0.L0); 
    solution_expand.E0 = LinearInterp(Par, Par_expand, solution0.E0); 
    
    solution_expand.L1 = LinearInterp(Par, Par_expand, solution0.L1); 
    solution_expand.E1 = LinearInterp(Par, Par_expand, solution0.E1); 
    solution_expand.V1 = LinearInterp(Par, Par_expand, solution0.V1); 
    solution_expand.EV1= LinearInterp(Par, Par_expand, solution0.EV1); 
    solution_expand.N1 = LinearInterp(Par, Par_expand, solution0.N1);       
    solution_expand.B1 = LinearInterp(Par, Par_expand, solution0.B1); 
    solution_expand.C1 = LinearInterp(Par, Par_expand, solution0.C1); 
    solution_expand.D1 = LinearInterp(Par, Par_expand, solution0.D1); 
    solution_expand.G1 = LinearInterp(Par, Par_expand, solution0.G1); 
    solution_expand.r1 = LinearInterp(Par, Par_expand, solution0.r1); 
    solution_expand.rd1= LinearInterp(Par, Par_expand, solution0.rd1);
    
    solution_expand.L1_Expand = LinearInterp(Par, Par_expand, solution0.L1_Expand); 
    solution_expand.E1_Expand = LinearInterp(Par, Par_expand, solution0.E1_Expand); 

    solution_expand.N1 = solution_expand.D1*(1-Par.theta) + solution_expand.E0 - solution_expand.L0 - solution_expand.B1;
    solution_expand.N1 = solution_expand.N1.*(solution_expand.N1<0);

  
    solution_expand.Profit_total   = LinearInterp(Par, Par_expand, solution0.Profit_total); 
    solution_expand.Profit_deposit = LinearInterp(Par, Par_expand, solution0.Profit_deposit); 
    solution_expand.Profit_loan    = LinearInterp(Par, Par_expand, solution0.Profit_loan); 
    solution_expand.Profit_finance = LinearInterp(Par, Par_expand, solution0.Profit_finance); 

    
    solution_expand.r_bar  = interp1(Par.fGrid,solution0.r_bar,Par_expand.fGrid);
    solution_expand.rd_bar = interp1(Par.fGrid,solution0.r_bar',Par_expand.fGrid); solution_expand.rd_bar= solution_expand.rd_bar';
    solution_expand.D_avg  = interp1(Par.fGrid,solution0.D_avg,Par_expand.fGrid);
    solution_expand.B_avg  = interp1(Par.fGrid,solution0.B_avg,Par_expand.fGrid);
    solution_expand.rd_star= interp1(Par.fGrid,solution0.rd_star,Par_expand.fGrid);
    solution_expand.D_star = interp1(Par.fGrid,solution0.D_star,Par_expand.fGrid);
    solution_expand.r_star = interp1(Par.fGrid,solution0.r_star,Par_expand.fGrid);
    solution_expand.B_star = interp1(Par.fGrid,solution0.B_star,Par_expand.fGrid);
    
    
    
    % Interpolate distribution ############################################
    
%     cum_distr = cumsum(reshape(solution0.distr,Par.nA*Par.nf,[]));
%     cum_distr_expand = LinearInterp(Par, Par_expand, cum_distr); 
%     cum_distr_expand = cumsum(reshape(cum_distr_expand,1,[]));
% 
%     solution_expand.distr = diff([0,cum_distr_expand]);
    
end
