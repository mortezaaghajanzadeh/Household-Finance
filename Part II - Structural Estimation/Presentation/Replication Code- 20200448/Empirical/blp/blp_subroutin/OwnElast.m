function [OwnSemiElast]=OwnElast(alpha,theta,delta,Data)

% ElastNestedLogit - Calculate price elasticities for Nested Logit models
% Syntax:  f=ElastLogit(alpha,theta,delta,Data)
%
% Inputs:
%    input1 - alpha     = price coefficient
%    input2 - thetaNL   = random coefficient
%    input1 - delta     = mean value
%    input2 - Data
%
% Outputs:
%    prods x prods elasticity matrix
%
% Subfunctions: NLShareCalculation 



% Unpack
qweight     = Data.qweight;
qv = Data.qv;
alpha_i = alpha+qv*theta(1);

% Market Shares
[sij,sj] = ShareCalculation(theta,delta,Data);

sij_alpha_i = sij.*alpha_i;
part1            = sij_alpha_i;          % diagonal in multiple dimensions
part3            = sij.*sij_alpha_i;

derSharerij=(part1 - part3);
derSharerj    = sum(derSharerij.*qweight,2);
OwnSemiElast = derSharerj./sj;



end
