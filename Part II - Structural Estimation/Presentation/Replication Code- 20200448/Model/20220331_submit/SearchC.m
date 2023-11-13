function [C1,E1,EV1,v] = SearchC(Par,E0,Profit,EV0,dummy_Kconst)

% This fucntion solves for the optimal dividends
% E0 is the initial book equity 
% L0 is the initial loan amount
% L1 is the end-of-period book equity

% Profit is the optimal profit (given firms' optimal choice of D)
% EV is the continuation value
% dummy_Kconst indicates whether the capital constraint is binding


if  size(E0) == size(dummy_Kconst) & size(E0) == size(Profit)
else 
    disp('matrix dimensions do not match in SearchD')
end

[n1,n2] = size(E0);

E0 = reshape(E0,1,[]);
Profit = reshape(Profit,1,[]);
dummy_Kconst = reshape(dummy_Kconst,1,[]);

E0_mat = repmat(E0,Par.nE_Expand,1);
Profit_mat= repmat(Profit,Par.nE_Expand,1);
dummy_mat = repmat(dummy_Kconst,Par.nE_Expand,1);

E1_mat = repmat(Par.EGrid_Expand',1,Par.nE*Par.nL*Par.nL_Expand);
EV_mat = reshape(EV0,[],Par.nE)';

if length(Par.EGrid_Expand) == length(Par.EGrid) & max(abs(Par.EGrid_Expand-Par.EGrid))<1e-4
else
    EV_mat = interp1(Par.EGrid,EV_mat,Par.EGrid_Expand,'linear','extrap');
end
EV_mat = repmat(EV_mat,1,Par.nE);

C_mat = E0_mat + Profit_mat *(1- Par.tauc) - E1_mat;
V_mat = C_mat + EV_mat;
V_mat = V_mat + (Par.C_low - C_mat).*(C_mat > Par.C_low).*dummy_mat;

V_mat = V_mat - Par.issue*(C_mat < Par.C_low);

[v,indexE1_Expand] = max(V_mat); 
E1 = Par.EGrid_Expand(indexE1_Expand);
C1 = E0 + Profit*(1 - Par.tauc) - E1; 
C1 = C1 + (Par.C_low - C1).*(C1 > Par.C_low).*dummy_Kconst;

index_mat  = ([1:n1*n2]-1)*Par.nE_Expand + indexE1_Expand;        
EV1 = v - C1;

v  = reshape(v,n1,n2); 
C1 = reshape(C1,n1,n2);
E1 = reshape(E1,n1,n2);
EV1= reshape(EV1,n1,n2);

end