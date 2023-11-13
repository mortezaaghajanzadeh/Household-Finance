function [r_star,B_star,Prof_loan] = calc_eqB_BLP(Par, f, Er)

    % If Er is not supplied
    % This function calculates the equilibrium loan rate by assuming that banks
    % are symmetric and equating bank's FOC on the loan market to 0
    
    % If Erd is supplied
    % This function calculates the optimal deposit rate by assuming that
    % other banks have an equivalent rate of Er

    r_Grid = min(f):0.00001:1;
    f      = reshape(f,1,[]); nf = length(f);
    r_Grid = reshape(r_Grid,[],1); nr = length(r_Grid);

    Ef = interp1(Par.fGrid, Par.Ef_Grid, f);
    
    f_Grid = repmat(f,nr,1);
    Ef_Grid= repmat(Ef,nr,1);
    r_Grid = repmat(r_Grid,1,nf);
    
    if nargin == 2       
        Er = r_Grid;        
    elseif nargin == 3       
        Er = reshape(kron(Er,ones(1,nr)),size(r_Grid));
    end
       
   
    exp_corporate  = exp(Par.nuC - Par.alpha *(Ef_Grid + 0*Par.meanA + (Ef_Grid-Par.meanf)*Par.rhoAf1));
    exp_0          = 1;
    exp_bank       = exp(Par.nuB - Par.alpha *Er)*(Par.NB-1);
    exp_banki      = exp(Par.nuB - Par.alpha *r_Grid);
    
    B_Grid = exp_banki./(exp_banki + exp_bank + exp_0 + exp_corporate);
    
%   PROFIT =  B*a0*¡¾r_mat - Ef¡¿
    FOC = Par.Eyrs*(r_Grid-Ef_Grid) - (Par.meanA+(Ef_Grid-Par.meanf)*Par.rhoAf1)*Par.Eyrs - Par.cb0*Par.Eyrs...
        - 1./(Par.alpha*(1-B_Grid))*Par.Eyrs;
    
    index0 = sum(FOC < 0); index0 = index0 + 1*(index0 == 0);
    index = ((1:nf)-1)*nr + index0;
    FOC0 = FOC(index);
    FOC1 = FOC(index+1);
    
    if FOC1.*FOC0 < 1E-4
    else
        disp('solution in calc_eqB_BLP is a corner solution')
    end

    r_star = r_Grid(index);
    B_star = B_Grid(index);
    
    Prof_loan = B_star.*(r_star - Ef)*Par.Eyrs - B_star.*(Par.meanA+(Ef-Par.meanf)*Par.rhoAf1)*Par.Eyrs - Par.cb0*B_star*Par.Eyrs;
end        

