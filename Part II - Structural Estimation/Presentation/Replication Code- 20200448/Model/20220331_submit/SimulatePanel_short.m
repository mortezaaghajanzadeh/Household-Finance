 function [firm, simu_con, simu_un, rec1] = SimulatePanel_short(Par, solution, shock)


% This function simulate a panel of banks using the solution and the shocks supplied 

   [t,N] = size(shock.valueA);
   
   indexfA0= zeros(t,N);
   index0  = zeros(t,N);
   
   indexE1 = zeros(t+1,N);
   indexL1 = zeros(t+1,N);
   indexEL1= zeros(t+1,N);
         
   indexE1(1,:) = shock.indexE(1,:);
   indexL1(1,:) = shock.indexE(1,:);
   indexEL1(1,:)= (indexE1(1,:) -1)*Par.nL+ indexL1(1,:);
      
   lgrid_temp = (Par.LGrid(1:end-1) + Par.LGrid(2:end))/2;
   egrid_temp = (Par.EGrid(1:end-1) + Par.EGrid(2:end))/2;

   temp_L1 = Par.LGrid(indexL1(1,:));
   temp_E1 = Par.EGrid(indexE1(1,:));
   
   for i = 1:t
       indexfA0(i,:) = (shock.indexf(i) -1)*Par.nA+ shock.indexA(i,:);
       if i > 1 
           
           temp_indexL1 = temp_L1*0+1;
           temp_indexE1 = temp_E1*0+1;
           
           for j = 1:Par.nL-1
              temp_indexL1 = temp_indexL1 + (temp_L1 > lgrid_temp(j));
           end
        
           for j = 1:Par.nE-1
              temp_indexE1 = temp_indexE1 + (temp_E1 > egrid_temp(j));
           end    
           
           indexE1(i,:) = temp_indexE1;    
           indexL1(i,:) = temp_indexL1; 
       end

          
       indexEL1(i,:)= (indexE1(i,:) -1)*Par.nL+ indexL1(i,:);
       index0(i,:) = (indexEL1(i,:)-1)*Par.nA*Par.nf + indexfA0(i,:);
                     
       temp_E1 = solution.E1_Expand(index0(i,:));
       temp_L1 = solution.L1_Expand(index0(i,:));
              
       exit_dummy = (solution.V1(index0(i,:)) < Par.V_low);
       temp_E1 = temp_E1 + (shock.valueE(i,:) - temp_E1).*exit_dummy;
       temp_L1 = temp_L1 + (shock.valueL(i,:) - temp_L1).*exit_dummy;
       
   end
   
   indexE0 = indexE1(1:t,:);
   indexL0 = indexL1(1:t,:);
   indexEL0= indexEL1(1:t,:);
   
   indexE1 = indexE1(2:t+1,:);
   indexL1 = indexL1(2:t+1,:);
   indexEL1= indexEL1(2:t+1,:);

   firm.A0 = solution.A0(index0);   
   firm.f0 = solution.f0(index0);
   firm.Ef0= solution.Ef0(index0);   
   firm.f00  = [firm.f0(1,:);firm.f0(1:end-1,:)];
   firm.f000 = [firm.f00(1,:);firm.f00(1:end-1,:)];
   firm.f0000= [firm.f000(1,:);firm.f000(1:end-1,:)];
   
   firm.df0  = repmat(shock.deltaf,1,N);
   firm.df00 = [firm.df0(1,:);firm.df0(1:end-1,:)];
   firm.df000= [firm.df00(1,:);firm.df00(1:end-1,:)];
   
   firm.L0 = solution.L0(index0);
   firm.L00 = [firm.L0(1,:);firm.L0(1:end-1,:)];
   firm.L000= [firm.L00(1,:);firm.L00(1:end-1,:)];
   
   firm.E0 = solution.E0(index0);
   if Par.mu > 0 firm.L1 = solution.L1_Expand(index0)/(1-Par.mu); else firm.L1 = solution.B1(index0); end
   firm.E1 = solution.E1_Expand(index0);
   firm.V1 = solution.V1(index0);
   firm.EV1= solution.EV1(index0);
   firm.N1 = solution.N1(index0);
   firm.B1 = solution.B1(index0);
   firm.C1 = solution.C1(index0);
   firm.D1 = solution.D1(index0);
   firm.G1 = solution.G1(index0);
   firm.r1 = solution.r1(index0);
   firm.rd1= solution.rd1(index0);
   
   firm. EV0 = [firm.EV1(1,:);firm.EV1(1:end-1,:)];
   
   firm.Profit_total   = solution.Profit_total(index0);
   firm.Profit_deposit = solution.Profit_deposit(index0);
   firm.Profit_loan    = solution.Profit_loan(index0);
   firm.Profit_finance = solution.Profit_finance(index0);   
      
   firm.indexE0 = indexE0;
   firm.indexL0 = indexL0;
   firm.indexA0 = shock.indexA;
   firm.indexf0 = repmat(shock.indexf,1,N);
   firm.indexEL0= indexEL0;
   firm.indexfA0= indexfA0;
   firm.index0  = index0;
   
   firm.indexE1 = indexE1;
   firm.indexL1 = indexL1;
   firm.indexEL1= indexEL1;
   
   
   
   % ############################ SIMULATE ################################
    firm.r_bar = solution.r_bar(firm.indexf0);
    firm.r2(1,:) = firm.r1(1,:);
    for i = 2:t
        firm.r2(i,:) = (firm.r2(i-1,:).*firm.L0(i,:) + firm.r1(i,:).*firm.B1(i,:))./(firm.L0(i,:)+firm.B1(i,:));
    end
    
    [~, SL_bank, SL_bond, SL_none] = calc_Bi_BLP(Par,firm.r1,firm.f0,firm.r_bar,'r'); 
    firm.T1 = 1-SL_none;
%     firm.T1(1,:) = SL_bond(1,:)/Par.mu + mean(firm.L1(i,:))*Par.NB; 
%     for i = 2:t
%         firm.T1(i,:) = firm.T1(i-1,:)*(1-Par.mu) + (SL_bond(i,:)+ mean(firm.L1(i,:))*Par.NB);
%     end
    
    f0 = firm.f0;
    f00= firm.f00;
    f000=firm.f000;
    f0000=firm.f0000;
    
    df0 = firm.df0;
    df00= firm.df00;
    df000=firm.df000;
    
    L0 = firm.L0;
    L00= firm.L00;
    L000=firm.L000;
    
    A0 = firm.A0; 
    B1 = firm.B1;
    L1 = firm.L1;
    E0 = firm.E0; 
    E1 = firm.E1; 
    V1 = firm.V1;
    T1 = firm.T1;
    N1 = firm.N1;
    G1 = firm.G1;
    C1 = firm.C1;
    D1 = firm.D1;
    EV1= firm.EV1;
    EV0 = firm.EV0;
    
    Profit_deposit = firm.Profit_deposit;
    Profit_loan = firm.Profit_loan;
    Profit_finance = firm.Profit_finance;
    Profit_total = firm.Profit_total;
    
    Ef0 = interp1(Par.fGrid, Par.Ef_Grid,firm.f0);
    rd1 = firm.rd1; 
    r1  = firm.r1;
    r2  = firm.r2;
    
    
    % Calculate at the aggregate level
    weights = A0./repmat(mean(A0,2),1,N);
    agg_B1 = mean(B1,2);
    agg_E1 = mean(E1,2); 
    agg_E0 = mean(E0,2); 
    agg_L1 = mean(L1,2); 
    agg_T1 = mean(T1,2); 
    agg_D1 = mean(D1,2);
    agg_N1 = mean(N1,2);
    agg_G1 = mean(G1,2);
    agg_A1 = mean(D1 - N1 + E1,2);
    
    agg_rd1= mean(weights.*rd1,2); 
    agg_r1 = mean(weights.*r1,2); 
    agg_r2 = mean(weights.*r2,2); 
    agg_f0 = mean(f0,2);
    agg_Ef0= mean(Ef0,2);
    agg_f00= mean(f00,2);
    agg_f000=mean(f000,2);
    agg_f0000=mean(f0000,2);
    
    agg_df0  =mean(df0,2);
    agg_df00 =mean(df00,2);
    agg_df000=mean(df000,2);
    
    agg_L0  =mean(L0,2);
    agg_L00 =mean(L00,2);
    agg_L000=mean(L000,2);


    a2 = regress(agg_B1,[agg_f0,agg_f0*0+1]); 
    b2 = regress(agg_L1,[agg_f0,agg_f0*0+1]); 
    c2 = regress(agg_T1,[agg_f0,agg_f0*0+1]);    
    d2 = regress(agg_f0-agg_rd1,[agg_f0,agg_f0*0+1]);    
    e2 = regress(agg_r1-agg_Ef0,[agg_f0,agg_f0*0+1]); 
    f2 = regress(agg_rd1,[agg_f0,agg_f0*0+1]);
    g2 = regress(agg_r2,[agg_f0,agg_f0*0+1]);
    h2 = regress(agg_D1,[agg_f0,agg_f0*0+1]);
    i2 = regress(agg_N1,[agg_f0,agg_f0*0+1]);
    j2 = regress(agg_E1,[agg_f0,agg_f0*0+1]);
    jj2 = regress(agg_E0,[agg_f0,agg_f0*0+1]);
    k2 = regress(agg_G1+agg_D1*Par.theta,[agg_f0,agg_f0*0+1]);
    l2 = regress(agg_A1,[agg_f0,agg_f0*0+1]);
    
    DSS_rate  = regress((agg_rd1(2:end)-agg_rd1(1:end-1)),[agg_f0(2:end)-agg_f0(1:end-1),agg_f0(1:end-1)*0+1]);
    DSS_spread   = regress((agg_f0(2:end)-agg_f0(1:end-1)) - (agg_rd1(2:end)-agg_rd1(1:end-1)),[agg_f0(2:end)-agg_f0(1:end-1),agg_f0(1:end-1)*0+1]);
    DSS_quantity = regress(log(agg_D1(2:end))-log(agg_D1(1:end-1)),[agg_f0(2:end)-agg_f0(1:end-1),agg_f0(1:end-1)*0+1]);
    
    SS_rate          = regress((agg_r1(2:end)-agg_r1(1:end-1)) ,[agg_f0(2:end)-agg_f0(1:end-1),agg_f0(1:end-1)*0+1]);
    SS_spread     = regress((agg_r1(2:end)-agg_r1(1:end-1))-(agg_Ef0(2:end)-agg_Ef0(1:end-1)) ,[agg_f0(2:end)-agg_f0(1:end-1),agg_f0(1:end-1)*0+1]);
    SS_quantity   = regress(log(agg_L1(2:end))-log(agg_L1(1:end-1)),[agg_f0(2:end)-agg_f0(1:end-1),agg_f0(1:end-1)*0+1]);

    simu_un.DSS_rate = DSS_rate(1);
    simu_un.DSS_spread = DSS_spread(1);
    simu_un.DSS_quantity =  DSS_quantity(1);
    simu_un.SS_rate = SS_rate(1);
    simu_un.SS_spread = SS_spread(1);
    simu_un.SS_quantity = SS_quantity(1);
    
    % Calculate at individul firm level
    dummy = (V1 < Par.V_low) + (B1 < 2*1E-4);
    
    f0(dummy >0) = [];  f0=reshape(f0,[],1);
    f00(dummy>0) = [];  f00=reshape(f00,[],1);
    f000(dummy>0) = []; f000=reshape(f000,[],1);
    f0000(dummy>0) = [];f0000=reshape(f0000,[],1);

    df0(dummy >0) = []; df0=reshape(df0,[],1);
    df00(dummy>0) = []; df00=reshape(df00,[],1);
    df000(dummy>0)= []; df000=reshape(df000,[],1);
    
    L0(dummy >0) = [];  L0=reshape(L0,[],1);
    L00(dummy>0) = [];  L00=reshape(L00,[],1);
    L000(dummy>0)= [];  L000=reshape(L000,[],1);
  
    B1(dummy >0) = [];  B1=reshape(B1,[],1);
    L1(dummy >0) = [];  L1=reshape(L1,[],1);
    E1(dummy >0) = [];  E1=reshape(E1,[],1);
    E0(dummy >0) = [];  E0=reshape(E0,[],1);
    A0(dummy >0) = [];  A0=reshape(A0,[],1);
    V1(dummy >0) = [];  V1=reshape(V1,[],1);
    T1(dummy >0) = [];  T1=reshape(T1,[],1);
    N1(dummy >0) = [];  N1=reshape(N1,[],1);
    G1(dummy >0) = [];  G1=reshape(G1,[],1);
    C1(dummy >0) = [];  C1=reshape(C1,[],1);
    D1(dummy >0) = [];  D1=reshape(D1,[],1);
    EV1(dummy >0) = []; EV1=reshape(EV1,[],1);
    EV0(dummy >0) = []; EV0=reshape(EV0,[],1); 
    
    r1(dummy >0) = [];  r1 =reshape(r1,[],1);
    r2(dummy >0) = [];  r2 =reshape(r2,[],1);
    rd1(dummy>0) = [];  rd1=reshape(rd1,[],1);
    Ef0(dummy>0) = [];  Ef0=reshape(Ef0,[],1);
   
    Profit_deposit(dummy>0) = []; Profit_deposit=reshape(Profit_deposit,[],1);
    Profit_loan(dummy>0) = []; Profit_loan = reshape(Profit_loan,[],1);
    Profit_finance(dummy>0) = []; Profit_finance = reshape(Profit_finance,[],1);
    Profit_total(dummy>0) = []; Profit_total = reshape(Profit_total,[],1);
     
    
     A1 = D1 - N1 + E1;
     EV1 = max(EV1,1e-4);
    
    simu_un.meanE = mean(E1);
    simu_un.meanL = mean(L1);
    simu_un.meanD = mean(D1);
    simu_un.meanB = mean(B1);
    simu_un.meanG = mean(G1);
    simu_un.meanV = mean(V1);
    simu_un.meanN = mean(N1);
    simu_un.meanC = mean(C1);
    simu_un.meanEV= mean(EV1);
    simu_un.meanT= mean(T1);
    
    simu_un.meanr   = mean(r1);
    simu_un.meanrd  = mean(rd1);
    simu_un.spreadrd= mean(f0-rd1);
    simu_un.spreadr = mean(-Ef0+r1-exp(Par.nuA+Par.sigA^2));
    
    simu_un.meanA = simu_un.meanD - simu_un.meanN + simu_un.meanE + simu_un.meanC;
    
    simu_un.meanL0 = mean(L0);
    simu_un.meanE0 = mean(E0);
    simu_un.meanA0 = mean(A0);
    simu_un.meanf0 = mean(f0);
    
    simu_un.meanPdeposit = mean(Profit_deposit);
    simu_un.meanPloan    = mean(Profit_loan);
    simu_un.meanPfinance = mean(Profit_finance);
    simu_un.meanPtotal   = mean(Profit_total);
    
    
    for i = 1:Par.nf
        indexi = find(firm.indexf0 == i);
        simu_con.r1(i)  = mean(firm.r1(indexi).*firm.B1(indexi))/mean(firm.B1(indexi));
        simu_con.rd1(i) = mean(firm.rd1(indexi).*firm.D1(indexi))/mean(firm.D1(indexi));
        simu_con.D1(i)  = mean(firm.D1(indexi));
        simu_con.B1(i)  = mean(firm.B1(indexi));
        simu_con.L1(i)  = mean(firm.L1(indexi));
        simu_con.L0(i)  = mean(firm.L0(indexi));
        simu_con.E0(i)  = mean(firm.E0(indexi));
        simu_con.N1(i)  = mean(firm.N1(indexi));
        simu_con.T1(i)  = mean(firm.T1(indexi));
        simu_con.profit(i)  = mean(firm.Profit_total(indexi));
    end

 
     simu_un.meanC2V = mean(A1.*C1./V1.*(C1./V1>0.0001).*(C1./V1<1))/mean(A1.*(C1./V1>0.001).*(C1./V1<1));
     simu_un.meanN2A = mean(A1.*N1./A1)/mean(A1);
     simu_un.meanN2D = mean(A1.*N1./D1)/mean(A1);
     simu_un.meanL2D = mean(A1.*L1./D1)/mean(A1);
     simu_un.meanD2A = mean(A1.*D1./A1)/mean(A1);
     simu_un.meanM2B = mean(A1.*min(max(V1./max(E1, 1e-4),1),10))/mean(A1);
     simu_un.meanLEV = mean(A1.*min(A1./max(E1,1e-4),20))/mean(A1);
     simu_un.meanFIX1 = mean(A1.*(D1*Par.cd0+(L1+B1)*Par.cb0 - (L1+B1)*Par.meanA+Par.fix+Par.a0.*abs(N1)+Par.a1*N1.^2./D1)./A1)/mean(A1);
     simu_un.meanFIX2 = mean(A1.*(D1*Par.cd0+(L1+B1)*Par.cb0 - (L1+B1)*Par.meanA+Par.fix)./A1)/mean(A1);
    
    
     simu_un.stdf0 = std(f0);
     simu_un.stdA0 = std(A0);
     simu_un.stdN2A  = (mean(A1.*(N1./A1).^2)/mean(A1)-simu_un.meanN2A^2)^0.5;
     simu_un.stdN2D  = (mean(A1.*(N1./D1).^2)/mean(A1)-simu_un.meanN2D^2)^0.5;

     a1 = regress(B1,[f0,f0*0+1]);
     b1 = regress(L1,[f0,f0*0+1]);      
     c1 = regress(T1,[f0,f0*0+1]);
     d1 = regress(f0-rd1,[f0,f0*0+1]);    
     e1 = regress(r1-Ef0,[f0,f0*0+1]);
     f1 = regress(rd1,[f0,f0*0+1]);    
     g1 = regress(r2,[f0,f0*0+1]); 
     h1 = regress(D1,[f0,f0*0+1]); 
     i1  = regress(N1, [f0,f0*0+1]);    
     j1  = regress(E1, [f0,f0*0+1]);
     jj1 = regress(E0, [f0,f0*0+1]);
     k1 = regress(G1+D1*Par.theta, [f0,f0*0+1]);
     l1 = regress(D1-N1+E1, [f0,f0*0+1]);


     % calculate dB/df (dT/df) scaled by total loans (credit) outstanding, aggregated over a two year horizon
     x = a1(1)/simu_un.meanB*Par.mu; rho = regress(agg_f0, agg_f00);  y =  x*rho + x*(1-Par.mu)^1;  
     simu_un.reg_B_agg_1y  = x;  
     simu_un.reg_B_agg  = y;  

     x = c1(1)/simu_un.meanT*Par.mu; rho = regress(agg_f0, agg_f00);  y =  x*rho + x*(1-Par.mu)^1;  
     simu_un.reg_T_agg_1y  = x;  
     simu_un.reg_T_agg  = y;  
     
     simu_un.reg_srd         = d1(1);
     simu_un.reg_srd_agg = d2(1);
     simu_un.reg_sr           = e1(1);
     simu_un.reg_sr_agg   = e2(1);
     simu_un.reg_rd          = f1(1);
     simu_un.reg_rd_agg  = f2(1);
     simu_un.reg_r            = g1(1);
     simu_un.reg_r_agg    = g2(1);

     
     simu_un.reg_L          = b1(1);
     simu_un.reg_L_agg  = b2(1);
     simu_un.reg_L2L      = b1(1)/simu_un.meanL;
     simu_un.reg_L2L_agg= b2(1)/simu_un.meanL;    
     simu_un.reg_D          = h1(1)/simu_un.meanD;
     simu_un.reg_D_agg  = h2(1)/simu_un.meanD;     
     simu_un.reg_N          = i1(1)/abs(simu_un.meanN);
     simu_un.reg_N_agg  = i2(1)/abs(simu_un.meanN);     
     simu_un.reg_E1         = j1(1)/abs(simu_un.meanE);
     simu_un.reg_E1_agg = j2(1)/abs(simu_un.meanE);     
     simu_un.reg_E0         = jj1(1)/abs(simu_un.meanE);
     simu_un.reg_E0_agg = jj2(1)/abs(simu_un.meanE);     
     simu_un.reg_G1         = k1(1)/abs(simu_un.meanG+simu_un.meanD*Par.theta);
     simu_un.reg_G1_agg = k2(1)/abs(simu_un.meanG+simu_un.meanD*Par.theta);     
     simu_un.reg_A1         = l1(1)/abs(simu_un.meanD-simu_un.meanN+simu_un.meanE);
     simu_un.reg_A1_agg = l2(1)/abs(simu_un.meanD-simu_un.meanN+simu_un.meanE);     
     
     temp = lscov([df0,f0*0+1],V1,A1);  car1 = temp(1)*0.01/simu_un.meanV;     
     simu_un.CAR1 = car1;
                          
     rec1 = [simu_un.spreadrd, simu_un.spreadr,simu_un.meanD/simu_un.meanA, simu_un.reg_B_agg, simu_un.reg_T_agg, simu_un.reg_D_agg,simu_un.reg_srd_agg,simu_un.reg_sr_agg];
     rec1 = full(rec1);

end

