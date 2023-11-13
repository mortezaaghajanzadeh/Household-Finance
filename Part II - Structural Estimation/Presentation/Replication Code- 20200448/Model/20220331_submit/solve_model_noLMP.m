function [solution] = solve_model_noLMP(Par)

             

         rd_bar0 = repmat(Par.rd_star,Par.nq,1);
         r_bar0  = Par.r_star;
         D_avg0  = Par.D_star;
         B_avg0  = Par.B_star;
         V0 = repmat(kron(ones(Par.nA,1),Par.Prof_loanable'*(1- Par.tauc))+kron(Par.Prof_deposit'*(1- Par.tauc),ones(Par.nA,1)),1,Par.nE*Par.nL)*Par.beta/(1-Par.beta) ...
                   + repmat(kron(Par.EGrid,ones(1,Par.nL)),Par.nA*Par.nf,1);

       
       
%         Par.L_high = max(Par.B_comp)/Par.mu*(1-Par.mu);
%         Par.L_low  = 1e-4;
%         Par.LExpand_low  = Par.L_low;
%         Par.LExpand_high = Par.L_high;
% 
%         Par.B_low  = 1e-4;
%         Par.B_high = max(Par.B_comp);
%         Par.E_low  = Par.L_high*0.01+max(Par.fix,0);
%         Par.E_high = Par.L_high*0.15;
% 
%         Par.EExpand_low  = Par.L_high*0.01;
%         Par.EExpand_high = Par.L_high*0.20;
% 
%         Par.C_low  = 0;  
%         Par.V_low  = -1;
% 
%         
%         Par.B_mid  = Par.B_low*0.75+Par.B_high*0.25;
%         Btemp1     = linspace(log(Par.B_low), log(Par.B_mid), (Par.nL_Expand+1)/2); Btemp1  = exp(Btemp1);
%         Btemp2     = linspace(Par.B_mid, Par.B_high,(Par.nL_Expand+1)/2); 
%         Par.BGrid_Expand  = [Btemp1,Btemp2(2:end)];
% 
%         Par.LGrid = linspace(Par.L_low,Par.L_high,Par.nL);
%         Par.LGrid_Expand = linspace(Par.LExpand_low,Par.LExpand_high,Par.nL_Expand);  
% 
%         Par.E_mid  = Par.E_low*0.75+Par.E_high*0.25;
%         Etemp1     = linspace(log(Par.E_low), log(Par.E_mid), (Par.nE+1)/2); Etemp1  = exp(Etemp1);
%         Etemp2     = linspace(Par.E_mid, Par.E_high,(Par.nE+1)/2); 
%         Par.EGrid  = [Etemp1,Etemp2(2:end)];
% 
% 
%         Par.EExpand_mid  = Par.EExpand_low*0.75+Par.EExpand_high*0.25;
%         Etemp1     = linspace(log(Par.EExpand_low), log(Par.EExpand_mid), (Par.nE_Expand+1)/2); Etemp1  = exp(Etemp1);
%         Etemp2     = linspace(Par.EExpand_mid, Par.EExpand_high,(Par.nE_Expand+1)/2); 
%         Par.EGrid_Expand  = [Etemp1,Etemp2(2:end)];

% ###########################################################################




    
    Dist_ELfA0 = kron(ones(1,Par.nE*Par.nL)/(Par.nE*Par.nL), ...
                      ones(1,Par.nf*Par.nA)/(Par.nf*Par.nA)*Par.fA_Trans^100);
 
    [rd_star0,D_star0] = calc_eqD_BLP(Par, Par.fGrid, rd_bar0);
    [r_star0, B_star0] = calc_eqB_BLP(Par, Par.fGrid, r_bar0);
  
    if max(abs(D_star0-Par.D_star)) < 0.01 && max(abs(rd_star0-Par.rd_star)) < 0.01
    else disp('Deposit Market EQ Solution is unstable')
        rd_star0 = Par.rd_star;
        D_star0 = Par.D_star;
    end
    
    Vdummy    = 0;
    Ddummy    = 0;
        
    EGrid = kron(Par.EGrid,ones(1,Par.nL));
    LGrid = kron(ones(1,Par.nE),Par.LGrid);
    fGrid = kron(Par.fGrid,ones(1,Par.nA)); EfGrid = kron(Par.Ef_Grid,ones(1,Par.nA));
    AGrid = kron(ones(1,Par.nf),Par.AGrid);
    
    fA_distribution = (ones(1,Par.nA*Par.nf)/(Par.nA*Par.nf))*Par.fA_Trans^1000;
    fA_distribution = reshape(fA_distribution,Par.nA,Par.nf);
    A_distribution  = reshape(sum(fA_distribution,2),1,[]);
    f_distribution  = reshape(sum(fA_distribution,1),1,[]);
    
    
    solution.nf  = Par.nf;
    solution.nA  = Par.nA;
    solution.nE  = Par.nE;
    solution.nL  = Par.nL;
    
    solution.fGrid  = Par.fGrid;
    solution.AGrid  = Par.AGrid;
    solution.EGrid  = Par.EGrid;
    solution.LGrid  = Par.LGrid;
    solution.Ef_Grid= Par.Ef_Grid;
    solution.EA_Grid= Par.EA_Grid;
    
    solution.fA_Trans= Par.fA_Trans;
    solution.A_Trans = Par.A_Trans;
    solution.f_Trans = Par.f_Trans;
    solution.A_distr = Par.A_distr;
    solution.f_distr = Par.f_distr;
    solution.meanA   = Par.meanA;
    solution.meanf   = Par.meanf;
   
    
    V1     = V0*0;
    C1     = V1;
    D1     = V1;
    rd1    = V1;
    B1     = V1;
    r1     = V1;
    isValid = V1;

    Profit_loan   = V1;
    Profit_total  = V1;
    Profit_deposit= V1;
    Profit_finance= V1;    
    Kconst = V1;
    Policy = V1;
    PolicyL= V1;
    ValueL_Expand= V1;
    ValueE_Expand= V1;
    PolicyE= V1;

    
    % RECORD VALUE FUNCTION ITERATIONS     
    record_V = [];
    record_D = [];
    record_C = [];
    record_Policy = [];
    record_distribution = [];
    record_rate = [];
    Kconst_mat_save = [];
        
    
    err_outer = 1;
    iter_outer= 0;
    while (max(err_outer) > 0.05)  &&  iter_outer < 5
        err_inner = 1;
        iter_inner= 0;
        
      % ADJUST THE PRECISION HERE
        while (err_inner > 0.01)
  
           EV0 = Par.fA_Trans*V0;
           for i = 1:Par.nf
               
               E0_mat= repmat(EGrid,Par.nL_Expand,1);
               L0_mat= repmat(LGrid,Par.nL_Expand,1);               
               B_mat = repmat(Par.BGrid_Expand',1,Par.nL*Par.nE);
                   
              [D_mat,rd_mat,B_mat,r_mat, G_mat, N_mat, Ptotal_mat,Ploan_mat,Pdeposit_mat,Pfinance_mat] ...
                   = SearchD_noLMP(Par,E0_mat,L0_mat,B_mat,i,1,r_bar0,rd_bar0,D_star0,rd_star0);             
             
               Ploan_rec_mat = Ploan_mat;
                   
               for j = 1:Par.nA                   
              
                  index0 = (i-1)*Par.nA + j;
                  
                  Ploan_mat =Par.beta*B_mat.*r_bar0(i)*Par.Eyrs - Par.beta*(L0_mat+B_mat)*Par.fGrid(i) ...
                                     - Par.beta*(L0_mat+ B_mat)*(Par.AGrid(j)+(Par.fGrid(i)-Par.meanf)*Par.rhoAf1)...
                                     - Par.beta*(L0_mat+ B_mat)*Par.cb0;                                                      
                  
                  Ptotal_mat = Ploan_mat + Pdeposit_mat + Pfinance_mat;  
                  
                  if Par.ksi > 0 && Par.mu < 1
                      load('Data\Kconst_save')
                      Kconst_mat = Kconst(index0,:);
                      Kconst_mat = repmat(Kconst_mat,Par.nL_Expand,1);
                 else 
                       Kconst_mat = B_mat*0 + 2;
                  end
                  

                  L1_mat= (L0_mat+ B_mat)*(1-Par.mu);
                
                  EV0i = Par.beta*(reshape(EV0(index0,:) ,Par.nL,Par.nE));
                  EV0_mat = interp1(Par.LGrid,EV0i,reshape(L1_mat(:,1:Par.nL),[],1),'linear','extrap'); EV0_mat = reshape(EV0_mat,Par.nL_Expand,[]);

                  [C_mat,E1_mat,EV_mat,v_mat] = SearchC(Par,E0_mat,Ptotal_mat,EV0_mat,(B_mat>Kconst_mat+1e-4));                                                    
                
                  v_mat = v_mat- 1e4*(B_mat>Kconst_mat+1e-4).*(B_mat>Par.B_low+1e-4);

                  [v,indexL1_Expand] = max(v_mat); 
                  index_mat  = ([1:(Par.nE*Par.nL)]-1)*Par.nL_Expand + indexL1_Expand;        
                  
                  L1_Expand = L1_mat(index_mat);
                  E1_Expand = E1_mat(index_mat);   
                  
                  LGrid_temp = (Par.LGrid(1:end-1) + Par.LGrid(2:end))/2;                  
                  indexL1 = 1;
                  for i_temp = 1:Par.nL-1
                      indexL1 = indexL1 + (L1_Expand > LGrid_temp(i_temp));
                  end           
                  
                  EGrid_temp = (Par.EGrid(1:end-1) + Par.EGrid(2:end))/2;                  
                  indexE1 = 1;
                  for i_temp = 1:Par.nE-1
                      indexE1 = indexE1 + (E1_Expand > EGrid_temp(i_temp));
                  end 
                  
                  isValid(index0,:)= (1-(v<Par.V_low));      
                  v = max(C_mat(index_mat) + EV_mat(index_mat),Par.V_low);
                  
                   V1(index0,:) = v;
                   EV1(index0,:)= EV_mat(index_mat);
                   C1(index0,:) = C_mat(index_mat); 
                   D1(index0,:) = D_mat(index_mat);
                   G1(index0,:) = G_mat(index_mat);
                   N1(index0,:) = N_mat(index_mat);

                   B1(index0,:) = B_mat(index_mat);
                   r1(index0,:) = r_mat(index_mat);
                   rd1(index0,:) = rd_mat(index_mat);
                   Profit_total(index0,:)   = Ptotal_mat(index_mat);
                   Profit_deposit(index0,:) = Pdeposit_mat(index_mat);
                   Profit_loan(index0,:)    = Ploan_mat(index_mat);
                   Profit_finance(index0,:) = Pfinance_mat(index_mat);
                   Kconst(index0,:)  = Kconst_mat(index_mat);
                    
                   PolicyE(index0,:) = indexE1; 
                   PolicyL(index0,:) = indexL1; 
                   ValueL_Expand(index0,:) = L1_Expand;     
                   ValueE_Expand(index0,:) = E1_Expand;
                   Policy(index0,:)  = (PolicyE(index0,:) -1)*Par.nL+ indexL1;  
                   
                   r1(isValid(index0,:)==0) = Par.r_star(i);
                   rd1(isValid(index0,:)==0)= Par.rd_star(i);
                   B1(isValid(index0,:)==0) = 1E-4;
                   D1(isValid(index0,:)==0) = D_star0(i);  

               end
           end       
           iter_inner = iter_inner+1;
           err_inner  = max(max(abs(V1-V0)./(abs(V0)+abs(V1)+0.01)));      

           str   = [' Viter = ', num2str(iter_inner), ' , error = ', num2str(err_inner)];
           disp(str);

         % RECORD VALUE FUNCTION ITERATIONS 
           record_V = [record_V;iter_inner,reshape(V1,1,[])];
           record_D = [record_D;iter_inner,reshape(D1,1,[])];
           record_C = [record_C;iter_inner,reshape(C1,1,[])];
           record_Policy = [record_Policy;iter_inner,reshape(Policy,1,[])];

           
        % INSERT McQueen Porteus Bounds
          if iter_inner > 10
             b_upper = max(max(V1-V0));
             b_lower = min(min(V1-V0));
             V1 = V1 + (b_upper+b_lower)/2;
          end

          
          V0 = V1;
          if (sum(err_inner) > 1e50 || iter_inner > 100)
              Vdummy = 1;
              break;           
          end
        end
        
        
      E1 = Par.EGrid(PolicyE);
      L1 = Par.LGrid(PolicyL);
      R1 = Par.theta*D1;

      solution.V1 = V0; 
      solution.EV1= EV1;
      solution.C1 = C1;
      solution.D1 = D1;
      solution.E1 = E1;
      solution.L1 = L1;
      solution.B1 = B1;
      solution.G1 = G1;
      solution.N1 = N1;
      solution.R1 = R1;
      solution.r1 = r1;
      solution.rd1= rd1;

      solution.L0 = repmat(LGrid,Par.nf*Par.nA,1);
      solution.E0 = repmat(EGrid,Par.nf*Par.nA,1);
      solution.A0 = repmat(AGrid',1,Par.nE*Par.nL,1);
      solution.f0 = repmat(fGrid',1,Par.nE*Par.nL,1);
      solution.Ef0 = repmat(EfGrid',1,Par.nE*Par.nL,1);
      
      solution.E1_Expand = ValueE_Expand;
      solution.L1_Expand = ValueL_Expand;
      
      solution.Profit_total   = Profit_total;
      solution.Profit_loan    = Profit_loan;
      solution.Profit_deposit = Profit_deposit;
      solution.Profit_finance = Profit_finance;
       
      solution.Policy  = Policy;
      solution.PolicyE = PolicyE;
      solution.PolicyL = PolicyL;
      
      solution.Vdummy = Vdummy;
      solution.isValid= isValid;
      solution.Kconst = Kconst;
      solution.distr  = Dist_ELfA0;
      solution.distr_EL = full(reshape(sum(reshape(solution.distr,Par.nA*Par.nf,[])),Par.nL,[]));
      
      solution.r_bar = r_bar0;
      solution.rd_bar= rd_bar0;
      solution.D_avg = D_avg0;
      solution.B_avg = B_avg0;
      
      [solution.rd_star,solution.D_star] = calc_eqD_BLP(Par, Par.fGrid, solution.rd_bar);
      [solution.r_star, solution.B_star] = calc_eqB_BLP(Par, Par.fGrid, solution.r_bar);

      
%         if iter_outer > 1 && iter_inner == 1
%            break;
%         end
                     
       % CALCULATE STEADY STATE DISTRIBUTION USING TRANSITION MATRIX 
%          Dist_ELfA1 = Dist_ELfA0 + 1*(sdystt(Par,Policy,isValid,Dist_ELfA0,Transition_fAEL) - Dist_ELfA0);
%          Dist_ELfA0 = Dist_ELfA1;
         
%         [rd_bar1, D_avg1, rd_avg1] = forecast(Par, rd1, D1, Dist_ELfA1, 'D');
%         [r_bar1,  B_avg1, r_avg1] = forecast(Par, r1,  B1, Dist_ELfA1, 'B');
         
%        CHECK MARGINAL DISTRIBUTION
%        Dist_EL = ones(1,Par.nf*Par.nA)*reshape(Dist_ELfA1,Par.nf*Par.nA,[]);
%        Dist_fA = ones(1,Par.nE*Par.nL)*(reshape(Dist_ELfA1,Par.nf*Par.nA,[]))';

%        Dist_A_ELf  = reshape(Dist_ELfA1,Par.nA,[]);
%        Dist_f_ELA  = reshape(Dist_A_ELf',Par.nf,[]);
%        Dist_f      = Dist_f_ELA*ones(Par.nA*Par.nL*Par.nE,1); Dist_f = reshape(Dist_f,1,[]);
                   

       
       [rd_bar1, r_bar1, D_avg1, B_avg1] = forecast_panel(Par, solution);
       
       out = r_bar0;
       for i = 1:10
             [out, share_bank, share_bond, share_none] = calc_Bi_BLP(Par, B_avg1, Par.fGrid, out ,'B');
       end
       r_bar1  = out;
       r_bar1  = max(r_bar1, Par.r_comp+1e-4);
       
       [rd_star1,D_star1] = calc_eqD_BLP(Par, Par.fGrid, rd_bar1);
       [r_star1, B_star1] = calc_eqB_BLP(Par, Par.fGrid, r_bar1);
             
       err1 = max(max(abs(rd_bar1-rd_bar0)./(abs(rd_bar1)+0.01)));
       err2 = max(abs(r_bar1-r_bar0)  ./(abs(r_bar1)+0.01));
       err3 = max(abs(D_avg1-D_avg0)  ./(abs(D_avg1)+0.01));
       err4 = max(abs(B_avg1-B_avg0)    ./(abs(B_avg1)+0.01));
       err_outer  = [err1, err2];
 
       record_distribution = [record_distribution;iter_outer,reshape(Dist_ELfA0,1,[])];
       record_rate = [record_rate;iter_outer,...
                      err1, err2, err3, err4, reshape(B_avg1,1,[]),reshape(D_avg1,1,[]),reshape(rd_star1,1,[]),reshape(r_star1,1,[])];
       
       str   = [' Diter = ', num2str(iter_outer), ' , error = ', num2str(err_outer)];
       disp(str);  
       disp(rd_star0);
    
       iter_outer= iter_outer+1;
       
       D_star0  = D_star1;
       B_star0  = B_star1;
       rd_star0 = rd_star1;
       r_star0  = r_star1;
       
       rd_bar0  = rd_bar1;
       r_bar0   = r_bar1;   
       B_avg0   = B_avg1;
       D_avg0   = D_avg1;  
       
    end
    
    save('Data\innerloop_iter','record_V')
    save('Data\outerloop_iter','record_distribution','record_rate')  
    
    if isempty(Kconst_mat_save) 
    else save('Data\Kconst_mat_save','Kconst_mat_save')
    end

    
%     TAKE AVERAGE OF THE LAST TWO ITERATIONS
%     RESOLVE THE VALUE FUNCTION TO MAKE THINGS CONSISTENT

%      [solution] = take_average(Par, V0, (r_bar0 + r_bar1)/2, (rd_bar0 + rd_bar1)/2);

end

