function [sij, sj]=ShareCalculation(theta,delta,Data)

% NLShareCalculation - Nested Logit share calculation
% Syntax:   [sij sijg sj sjg s0 numer1 denom1 numer2 denom2]=NLShareCalculation(thetaNL,delta,Data)
%
% Inputs:
%    input1 - Non linear parameters
%    input2 - Data
%
% Outputs:
%    sij    = individual market shares
%    sijg   = individual conditional share of choosing good j from group g
%    sj     = market shares
%    sjg    = conditional share of choosing good j from group g 
%    s0     = outside good
%    numer1 = numerator1
%    denom1 = denominator1
%    numer2 = numerator 2
%    denom2 = denominator 2
%
% Subfunctions: none


% Unpack
qweight          = Data.qweight;
xv             = Data.xv;
nodes = Data.nodes;
nrc = Data.nrc;
dum_t = Data.dum_t;
id_t = Data.id_t;


% Market share calculation
mu=zeros(size(xv,1),nodes);
for i=1:nrc
    mu=mu+xv(:,(i-1)*nodes+1:i*nodes).*theta(i);
end
mudel           = bsxfun(@plus,delta,mu);
numer1          = exp(mudel);

% use multidimensional arrays
denom1=dum_t'*numer1+1;
denom1=denom1(id_t,:);

sij              = numer1./denom1;
sj                = sum(qweight.*sij,2);

end