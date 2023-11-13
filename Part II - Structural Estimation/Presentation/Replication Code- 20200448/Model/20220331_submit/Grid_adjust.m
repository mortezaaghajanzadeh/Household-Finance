function [Par] = Grid_adjust(Par)

% Adjust the range for FFR to deal with different FFR processes (in different subsamples)
% DEFINE THE GRIDS FOR DEPOSITOR HETEROGENEITY

    Par.qv = [0:.1:1]';
    Par.nq = length(Par.qv);
    Par.qweight = repmat(1/Par.nq,Par.nq,1);

%  [Par.qv, Par.qweight]=nwspgr('KPN', 1, 1);
%   Par.nq = length(Par.qweight);
%   Par.qv = exp(Par.qv);

% DEFINE THE GRIDS FOR SHOCK SPACE

    Par.f_low  = log(0.001);
    Par.f_high = Par.nuf+2*Par.sigf/(1-Par.rhof^2)^0.5;

    Par.f_mid  = log(exp(Par.f_high)/4);
    
    Par.A_low  = log(0.001);
    Par.A_high = log(0.030);
    Par.A_mid  = Par.A_low*0.75 + Par.A_high*0.25;

    Atemp1       = linspace(Par.A_low, Par.A_mid, (Par.nA+1)/2); Atemp1  = exp(Atemp1 );
    Atemp2       = linspace(exp(Par.A_mid), exp(Par.A_high),(Par.nA+1)/2); 
    
    Par.AGrid    = linspace(exp(Par.A_low), exp(Par.A_high), Par.nA); 
%     Par.AGrid    =  [Atemp1,Atemp2(2:end)];  if Par.nA ==1 Par.AGrid = exp(Par.nuA); end
                   
    ftemp1       = linspace(Par.f_low, Par.f_mid, (Par.nf+1)/2); ftemp1  = exp(ftemp1 );
    ftemp2       = linspace(exp(Par.f_mid), exp(Par.f_high),(Par.nf+1)/2); 
    Par.fGrid    = [ftemp1,ftemp2(2:end)];
    
%     Par.AGrid    = linspace(exp(Par.A_low), exp(Par.A_high), Par.nA); if Par.nA ==1 Par.AGrid = exp(Par.nuA); end
%     Par.fGrid    = linspace(exp(Par.f_low), exp(Par.f_high), Par.nf); 
    
    [Par.fA_Trans, Par.A_Trans, Par.f_Trans] = fATransition(Par);
%     [Par.fA_Trans, Par.A_Trans, Par.f_Trans] = CIR_Transition(Par);
    Par.ELfA = kron(ones(1,Par.nE*Par.nL)/(Par.nE*Par.nL),ones(1,Par.nf*Par.nA)/(Par.nf*Par.nA)*(Par.fA_Trans)^1000);
    
    Par.A_distr = ones(1,Par.nA)/(Par.nA)*(Par.A_Trans)^1000;
    Par.f_distr = ones(1,Par.nf)/(Par.nf)*(Par.f_Trans)^1000;
    Par.meanf = sum(Par.fGrid.*Par.f_distr);
    Par.meanA = sum(Par.AGrid.*Par.A_distr);
    
    

% DEFINE LONG-TERM RATES

  Ef = Par.fGrid';
  for i = 1:100
      Ef  =  Ef + (1 - Par.mu)^i*(Par.beta^i)*(Par.f_Trans^i)*Par.fGrid';
  end

  
  EA = Par.AGrid';
  for i = 1:100
      EA  =  EA + (1 - Par.mu)^i*(Par.beta^i)*(Par.A_Trans^i)*Par.AGrid';
  end
 
  Eyrs = 1;
  for i = 1:100
      Eyrs  =  Eyrs + (1 - Par.mu)^i*(Par.beta^i)*1;
  end
  
  Par.Ef_Grid = reshape(Ef/Eyrs,1,[]);
  Par.EA_Grid = reshape(EA/Eyrs,1,[]);
  Par.Eyrs = Eyrs;
  
  
Ef_simple = Par.fGrid';
  for i = 1:100
      Ef_simple  =  Ef_simple + (1 - Par.mu)^i*(Par.f_Trans^i)*Par.fGrid';
  end

  Eyrs_simple = 1;
  for i = 1:100
      Eyrs_simple  =  Eyrs_simple + (1 - Par.mu)^i*1;
  end
  
  Par.Ef_Grid_simple = reshape(Ef_simple/Eyrs_simple,1,[]);

  
  
% CALCULATE STEADY STATE

    Par.rd_low = 1E-4;
    Par.r_low  = 1E-4;

    [Par.rd_star,Par.D_star,Par.Prof_deposit] = calc_eqD_BLP(Par, Par.fGrid);
    [Par.r_star,Par.B_star,Par.Prof_loanable] = calc_eqB_BLP(Par, Par.fGrid);

     Par.r_comp = Par.Ef_Grid +(Par.meanA+(Par.Ef_Grid-Par.meanf)*Par.rhoAf1) + Par.cb0;
     Par.r_comp = Par.r_comp+0.007;
     Par.rd_comp= Par.fGrid - (Par.fGrid*Par.theta + Par.cd0);
     Par.D_comp = calc_Di_BLP_heterogenous(Par, Par.fGrid, Par.rd_comp, repmat(Par.rd_comp,Par.nq,1));
     Par.B_comp = calc_Bi_BLP(Par, Par.r_comp, Par.fGrid, Par.r_comp, 'r'); 
         
     
% Par.LGrid IS THE STATE SPACE. CONSTRUCTED AS [L-low, L-low/(1-Par.mu), L-low/(1-Par.mu)^2, ...]
% Par.LGrid_Expand IS THE CHOICE SPACE IF SEARCHING ON L, CONSTRUCTED AS: [LExpand-low, LExpand-low/(1-Par.mu), L-low/(1-Par.mu)^2, ...]
% Par.BGrid_Expand IS THE CHOICE SPACE IF SEARHCING ON B, CONSTRUCTED AS: linspace(B-low, B-high)


    Par.L_high = max(Par.B_star)/Par.mu*(1-Par.mu);
    Par.L_low  = 1e-4;
    Par.LExpand_low  = Par.L_low;
    Par.LExpand_high = Par.L_high;

    Par.B_low  = 1e-4;
    Par.B_high = max(Par.B_star);
    Par.E_low  = Par.L_high*0.01+max(Par.fix,0);
    Par.E_high = Par.L_high*0.15;
 
    Par.EExpand_low  = Par.L_high*0.01;
    Par.EExpand_high = Par.L_high*0.20;

    Par.C_low  = 0;  
    Par.V_low  = -1;

% ###########################################################################

    Par.B_mid  = Par.B_low*0.75+Par.B_high*0.25;
    Btemp1     = linspace(log(Par.B_low), log(Par.B_mid), (Par.nL_Expand+1)/2); Btemp1  = exp(Btemp1);
    Btemp2     = linspace(Par.B_mid, Par.B_high,(Par.nL_Expand+1)/2); 
    Par.BGrid_Expand  = [Btemp1,Btemp2(2:end)];
    
    Par.LGrid = linspace(Par.L_low,Par.L_high,Par.nL);
    Par.LGrid_Expand = linspace(Par.LExpand_low,Par.LExpand_high,Par.nL_Expand);  
    
    Par.E_mid  = Par.E_low*0.75+Par.E_high*0.25;
    Etemp1     = linspace(log(Par.E_low), log(Par.E_mid), (Par.nE+1)/2); Etemp1  = exp(Etemp1);
    Etemp2     = linspace(Par.E_mid, Par.E_high,(Par.nE+1)/2); 
    Par.EGrid  = [Etemp1,Etemp2(2:end)];


    Par.EExpand_mid  = Par.EExpand_low*0.75+Par.EExpand_high*0.25;
    Etemp1     = linspace(log(Par.EExpand_low), log(Par.EExpand_mid), (Par.nE_Expand+1)/2); Etemp1  = exp(Etemp1);
    Etemp2     = linspace(Par.EExpand_mid, Par.EExpand_high,(Par.nE_Expand+1)/2); 
    Par.EGrid_Expand  = [Etemp1,Etemp2(2:end)];

    
    
% ###########################################################################
     
%     Par.BGrid_Expand = linspace(Par.B_low, Par.B_high, Par.nL_Expand);  
%     Par.EGrid  = exp(linspace(log(Par.E_low), log(Par.E_high), Par.nE));
%     Par.EGrid_Expand  = exp(linspace(log(Par.EExpand_low), log(Par.EExpand_high), Par.nE_Expand));
%      

%      indexB = 1;
%      mB =(Par.BGrid_Expand(1:end-1)+Par.BGrid_Expand(2:end))/2;
%      for i = 1:length(Par.BGrid_Expand)-1
%          indexB = indexB + (Par.B_star > mB(i));
%      end
%      Par.B_star = Par.BGrid_Expand(indexB);
%      
%      indexB = 1;
%      mB =(Par.BGrid_Expand(1:end-1)+Par.BGrid_Expand(2:end))/2;
%      for i = 1:length(Par.BGrid_Expand)-1
%          indexB = indexB + (Par.B_comp > mB(i));
%      end
%      Par.B_comp =Par.B_comp + (Par.BGrid_Expand(indexB) - Par.B_comp).*(Par.B_comp<=Par.BGrid_Expand(end));

end