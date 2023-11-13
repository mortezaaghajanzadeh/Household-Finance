function dnd = diagnd(x)

% dnd = diagonal of multidimensional matrix
% Syntax:   dnd = diagnd(x)
%
% Inputs:
%   input1 - multidimensional matrix
%
% Outputs:
%   dnd = diagonal
%   Subfunctions: none
%

% Author: Laura Grigolon
% August 2012;

z   = size(x);
m   = z(1);
n   = z(2);
mn  = m * n;
p   = prod(z)/mn;
k   = (1:m+1:m*n)';
k   = repmat(k,1,p);
k   = bsxfun(@plus,k,(0:(p-1))*mn);
dnd = x(k);
end