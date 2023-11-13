 function [share_banki, banki_mk, share_cash, share_bond, share_bank] = calc_Di_BLP_heterogenous(Par,f,inp,Erd)

% This fucntion calcualtes the deposit shares given the deposit rates that banks offer

% f is the federal funds rate
% inp is the deposit rate offered by bank i, both size(f)==size(inp)
% Erd is the estimated deposit rate offered by other banks. Erd depends on
% the deteogenous deposit sensitivity of depositors size(Erd) =[Par.nq, length(f)]


    if size(inp) == size(f)
    else
            disp('matrix dimensions do not match in calc_Di_BLP_heterogenous')
    end
    
    if size(Erd) == [Par.nq, length(f)]
    else
            disp('matrix dimensions do not match in calc_Di_BLP_heterogenous')
    end

    mean0 = sum(Par.qv.*Par.qweight);
    sensitivity0 = Par.alphaD + Par.sigalphaD*(Par.qv - mean0); sensitivity0 = sensitivity0 + (1-sensitivity0).*(sensitivity0<1);

    nf = length(f);
    
    sensitivity0 = reshape(sensitivity0,[],1);
    density0     = reshape(Par.qweight,[],1);

    f   = reshape(f,1,[]);    
    rdi = reshape(inp,1,[]);   
    
    Erd = reshape(Erd,[],nf);   

    sensitivity_Grid = repmat(sensitivity0,1,nf);
    density_Grid = repmat(density0,1,nf);
    f_Grid   = repmat(f,  Par.nq,1);
    rdi_Grid = repmat(rdi,Par.nq,1);    
           
    
    % NORMALIZE EVERYTHING BY THE DEMAND OF MONEY
    log_cash  = sensitivity_Grid * 0 + Par.lM + Par.e;
    log_bond  = sensitivity_Grid .* f_Grid + Par.lB + Par.e;
    log_banki = sensitivity_Grid .* rdi_Grid + Par.lD + Par.e;
    log_Ebank = sensitivity_Grid .* Erd + Par.lD + 0;
    
    exp_cash  = exp(log_cash);
    exp_bond  = exp(log_bond);
    exp_banki = exp(log_banki);
    exp_Ebank = exp(log_Ebank); 
    
    
    if length(exp_cash(:,1)) == 1
        
        share_banki0 = (exp_banki./(exp_cash + (Par.ND-1)*exp_Ebank + exp_banki + exp_bond));
        share_cash0  = (exp_cash ./(exp_cash + (Par.ND-1)*exp_Ebank + exp_banki + exp_bond));
        share_bond0  = (exp_bond ./(exp_cash + (Par.ND-1)*exp_Ebank + exp_banki + exp_bond));
        share_Ebank0 = (exp_Ebank./(exp_cash + (Par.ND-1)*exp_Ebank + exp_banki + exp_bond))*(Par.ND-1); share_bank0 = share_Ebank0 + share_banki0;
  
        banki_FOC0 = sensitivity_Grid.*(exp_banki./(exp_cash + (Par.ND-1)*exp_Ebank + exp_banki + exp_bond)) ...
                                    .*(1-exp_banki./(exp_cash + (Par.ND-1)*exp_Ebank + exp_banki + exp_bond));

    else
    
        bankij = exp_banki./(exp_cash + (Par.ND-1)*exp_Ebank + exp_banki + exp_bond);
        cashj  = exp_cash ./(exp_cash + (Par.ND-1)*exp_Ebank + exp_banki + exp_bond);
        bondj  = exp_bond ./(exp_cash + (Par.ND-1)*exp_Ebank + exp_banki + exp_bond);
        Ebankj = exp_Ebank./(exp_cash + (Par.ND-1)*exp_Ebank + exp_banki + exp_bond)*(Par.ND-1);
      
        share_banki0 = sum(bankij.*density_Grid);
        share_cash0  = sum(cashj .*density_Grid);
        share_bond0  = sum(bondj .*density_Grid);
        share_Ebank0 = sum(Ebankj.*density_Grid); share_bank0 = share_Ebank0 + share_banki0;
    
        banki_FOC0 = sum(sensitivity_Grid.*density_Grid .*bankij .*(1-bankij));

    end
    
    
    
    log_cash  = sensitivity_Grid * 0 + Par.lM + Par.e;
    log_bond  = sensitivity_Grid .* f_Grid + Par.lB + Par.e;
    log_banki = sensitivity_Grid .* rdi_Grid + Par.lD ;
    log_Ebank = sensitivity_Grid .* Erd + Par.lD;
    log_Ebanke = sensitivity_Grid .* Erd + Par.lD + Par.e;

    exp_cash  = exp(log_cash);
    exp_bond  = exp(log_bond);
    exp_banki = exp(log_banki);
    if Par.ND > 1
        exp_Ebank = (exp(log_Ebank)*(Par.ND-2)+ exp(log_Ebanke))/(Par.ND-1);
    else
        exp_Ebank = 0;
    end

        
    if length(exp_cash(:,1)) == 1

        share_banki1 = (exp_banki./(exp_cash + (Par.ND-1)*exp_Ebank + exp_banki + exp_bond));
        share_cash1  = (exp_cash ./(exp_cash + (Par.ND-1)*exp_Ebank + exp_banki + exp_bond));
        share_bond1  = (exp_bond ./(exp_cash + (Par.ND-1)*exp_Ebank + exp_banki + exp_bond));
        share_Ebank1 = (exp_Ebank./(exp_cash + (Par.ND-1)*exp_Ebank + exp_banki + exp_bond))*(Par.ND-1); share_bank1 = share_Ebank1 + share_banki1;

        banki_FOC1 = sensitivity_Grid.*(exp_banki./(exp_cash + (Par.ND-1)*exp_Ebank + exp_banki + exp_bond)) ...
                                        .*(1-exp_banki./(exp_cash + (Par.ND-1)*exp_Ebank + exp_banki + exp_bond));

    else

        bankij = exp_banki./(exp_cash + (Par.ND-1)*exp_Ebank + exp_banki + exp_bond);
        cashj  = exp_cash ./(exp_cash + (Par.ND-1)*exp_Ebank + exp_banki + exp_bond);
        bondj  = exp_bond ./(exp_cash + (Par.ND-1)*exp_Ebank + exp_banki + exp_bond);
        Ebankj = exp_Ebank./(exp_cash + (Par.ND-1)*exp_Ebank + exp_banki + exp_bond)*(Par.ND-1);

        share_banki1 = sum(bankij.*density_Grid);
        share_cash1  = sum(cashj .*density_Grid);
        share_bond1  = sum(bondj .*density_Grid);
        share_Ebank1 = sum(Ebankj.*density_Grid); share_bank1 = share_Ebank1 + share_banki1;

        banki_FOC1 = sum(sensitivity_Grid.*density_Grid .*bankij .*(1-bankij));
            
    end

        share_banki  = (share_banki0 + share_banki1*(Par.ND-1))/Par.ND;
        share_cash   = (share_cash0  + share_cash1*(Par.ND-1)) /Par.ND;
        share_bond   = (share_bond0  + share_bond1*(Par.ND-1)) /Par.ND;
        share_bank   = (share_bank0  + share_bank1*(Par.ND-1)) /Par.ND;
        
        banki_mk     = share_banki./((banki_FOC0 + banki_FOC1*(Par.ND-1))/Par.ND);
    
end

