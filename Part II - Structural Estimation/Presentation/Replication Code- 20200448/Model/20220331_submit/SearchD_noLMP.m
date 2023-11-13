function [D1, rd1, B1, r1, G1, N1, profit0, profit1, profit2, profit3] = SearchD_noLMP(Par,E0,L0,B1,i,j,r_bar,rd_bar,D_star,rd_star)

% This fucntion solves for the optimal deposit rate and quantity
% E0 is the initial book equity 
% L0 is the initial loan amount
% L1 is the end-of-period book equity

% i and j are indeices for the demand shifter and the fed funds rate
% s_bar0 is the anticipated aggregate spread given f
% D_bar0 is the anticipated aggregate spread given f

% D1 is the optimal deposit quantity
% s1 is the optimal deposit spread

if size(E0) == size(B1) & size(E0) == size(L0)
else 
    disp('matrix dimensions do not match in SearchD')
end

[n1,n2] = size(E0);

f = Par.fGrid(i);
Ef = Par.Ef_Grid(i);
Er = r_bar(i);
Erd = rd_bar(:,i);

% Yifei note: The input values have sizes of (n_l by nl*n_e). The first
% dimension is future choice of l and the second dimension is the current
% state pair (e,l). The following code first rearrange all these into one
% long row vector. This row vector is then expanded into a matrix, with
% each row representing a choice of deposit. The row dimension is then
% maximized out.
E0 = reshape(E0,1,[]);
L0 = reshape(L0,1,[]);
B1 = reshape(B1,1,[]); B1(B1<0)=nan;
r1 = calc_Bi_BLP(Par,B1,f,Er,'B');  r1(r1<0)=nan;

% Amount lent out
P1 = B1*Par.mu./(Par.mu+r1);
P1 = B1;

% THERE IS NO CLOSED-FORM SOLUTION, USING INTERPOLATION TO CALCULATE RD(D)
n_steps = 9999;
%rd_grid = linspace(0.001,f,n_steps);
rd_grid = linspace(Par.rd_low,f*1.5,n_steps);
[D_grid, Dmk_grid]  = calc_Di_BLP_heterogenous(Par, repmat(f,1,n_steps), rd_grid, repmat(Erd,1,n_steps)); D_grid = Par.W0*D_grid;


% Expand permissible deposit taking amounts
D_upper = (L0+P1-E0)/(1-Par.theta);
D_upper = min(D_upper,min((2*Par.W0/Par.ND),1));
D_lower = max(D_star(i)*0.9,min(D_grid+1e-4));
D_step  = (D_upper-D_lower)/(Par.nD-1);

D_mat = [zeros(size(E0));kron([1:Par.nD-1]',D_step)];
D_mat = D_mat + D_lower;
D_mat = max(D_lower,D_mat)*0 + Par.D_star(i); rd_mat = D_mat*0 + Par.rd_star(i);




% [~, index_unique] = unique(D_grid);
% rd_mat  = interp1(D_grid(index_unique), rd_grid(index_unique), D_mat);
% Dmk_mat = interp1(D_grid(index_unique), Dmk_grid(index_unique), D_mat);

GN_mat = (1-Par.theta)*(1-Par.theta2)*D_mat - repmat( - E0 + P1 + L0, Par.nD, 1);
G_mat  = GN_mat.*(GN_mat>0) + (1-Par.theta2)*D_mat;
N_mat  = GN_mat.*(GN_mat<0);
L0_mat = repmat(L0, Par.nD,1);
r_mat  = repmat(r1, Par.nD,1);
B_mat  = repmat(B1, Par.nD,1);
P_mat  = repmat(P1, Par.nD,1);

temp_matrix = ((1-Par.theta)*D_mat*f - D_mat.*rd_mat - Par.cd0*D_mat - Par.cd1*D_mat.^2/2);
 [v_mat,index1] = max(temp_matrix);
    
profit_mat1 = P_mat.*(r_mat - Ef)*Par.Eyrs;

profit_mat2 = (1-Par.theta)*D_star(i)*f - D_star(i).*rd_star(i) - Par.cd0*D_star(i) - Par.cd1*D_star(i).^2/2 + repmat(E0*f,Par.nD,1) - Par.fix; 

profit_mat3 = ((1-Par.theta)*D_mat*f - D_mat.*rd_mat - Par.cd0*D_mat - Par.cd1*D_mat.^2/2) ...
             -((1-Par.theta)*D_star(i)*f - D_star(i).*rd_star(i) - Par.cd0*D_star(i) - Par.cd1*D_star(i).^2/2) ...
             + Par.a0 *N_mat - Par.a1*N_mat.^2/0.04/2; 
         
         
profit_mat1 = P_mat.*r_mat*Par.Eyrs - (L0_mat+P_mat).*f;
profit_mat2 = ((1-Par.theta)*D_mat*f - D_mat.*rd_mat - Par.cd0*D_mat - Par.cd1*D_mat.^2/2) + repmat(E0*f,Par.nD,1) - Par.fix; 
profit_mat3 =  Par.a0 *N_mat - Par.a1*N_mat.^2/Par.D_star(i)/2; 


profit_mat1 = profit_mat1*Par.beta;         
profit_mat2 = profit_mat2*Par.beta;         
profit_mat3 = profit_mat3*Par.beta;         
         
profit_mat0 = (profit_mat1 + profit_mat2 + profit_mat3);

        
[~,index1] = max(profit_mat0);
index_mat1 = ((1:n1*n2)-1)*Par.nD+index1;


% Using FOC
% FOC_deposit  =f*(1-Par.theta)- rd_mat - Par.cd0 - Par.cd1*D_mat - Dmk_mat;
% FOC_finance = (Par.a0  - Par.a1*GN_mat)*(1-Par.theta);
% FOC_finance =  FOC_finance.*(GN_mat < 0);
% FOC1   = abs(FOC_deposit + FOC_finance);
% [~,index1] = min(FOC1);
% index_mat1 = ((1:n1*n2)-1)*Par.nD+index1;



D1  = D_mat(index_mat1);
rd1 = rd_mat(index_mat1); 
G1  = G_mat(index_mat1);
N1  = N_mat(index_mat1);

profit0 = profit_mat0(index_mat1);
profit1 = profit_mat1(index_mat1);
profit2 = profit_mat2(index_mat1);
profit3 = profit_mat3(index_mat1);

G1  = reshape(G1,n1,n2);
N1  = reshape(N1,n1,n2);
D1  = reshape(D1,n1,n2);
rd1 = reshape(rd1,n1,n2);
B1  = reshape(B1,n1,n2);
r1  = reshape(r1,n1,n2);

profit0  = reshape(profit0, n1, n2);
profit1  = reshape(profit1, n1, n2);
profit2  = reshape(profit2, n1, n2);
profit3  = reshape(profit3, n1, n2);
end