function [rd_star,D_star,Prof_deposit] = calc_eqD_BLP(Par, f, Erd)
    
    % If Erd is not supplied
    % This function calculates the equilibrium deposit rate by assuming that
    % banks are symmetric and equating  abank's FOC on the deposit market to 0
    
    % If Erd is supplied
    % This function calculates the optimal deposit rate by assuming that
    % other banks have an equivalent rate of Erd 
    
    rd_Grid = Par.rd_low:1e-4:max(f);
    f  = reshape(f,1,[]); nf = length(f);
    rd_Grid = reshape(rd_Grid,[],1); nrd = length(rd_Grid);
    
    f_Grid  = repmat(f,nrd,1);
    rd_Grid = repmat(rd_Grid,1,nf);
    rd_Grid = min(rd_Grid,f_Grid);
    rd_Grid = max(rd_Grid,f_Grid-0.05);
        
    if nargin == 2       
        Erd = repmat(reshape(rd_Grid,1,[]),Par.nq,1);        
    elseif nargin == 3       
        Erd = kron(Erd,ones(1,nrd));
    end
    
    
    [D_Grid, D_mk] = calc_Di_BLP_heterogenous(Par,reshape(f_Grid,1,[]),reshape(rd_Grid,1,[]),Erd);
    D_Grid = reshape(D_Grid,nrd,nf);
    D_mk = reshape(D_mk,nrd,nf);

%    FOC = (f_Grid*(1-Par.theta)-rd_Grid-Par.cd0-Par.cd1*D_Grid).*dD_drd_Grid - D_Grid;
    FOC = (f_Grid*(1-Par.theta)-rd_Grid-Par.cd0-Par.cd1*D_Grid) - D_mk;
    
    [FOC0,index0] = min(abs(FOC)); index = ((1:nf)-1)*nrd + index0;
    
    rd_star = rd_Grid(index);
    D_star  = D_Grid(index)*Par.W0;
    
    if max(max(abs(FOC0))) > 1e-4 & rd_star ~= rd_Grid(1)
        disp('solution in calc_eqD_BLP does not converge')
    end

    Prof_deposit = D_star.* (f*(1-Par.theta) - rd_star - Par.cd0 - Par.cd1*D_star/2);
  
end 