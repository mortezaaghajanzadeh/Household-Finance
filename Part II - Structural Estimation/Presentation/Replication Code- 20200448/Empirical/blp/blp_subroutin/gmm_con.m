function [f , df]=gmm_con(theta,Data)
% gmmNL - This function is used to calculate the GMM estimator for the
% Nested Logit function
% 
% [f df]=gmmNL(thetaNL,BLPdata)
%
% Inputs:
%    input1 - Non linear parameters nested logit
%    input2 - Data
%
% Outputs:
%    f = gmm function value; df = gradient of f 
%
% Subfunctions: deltaLogit(theta,BLPdata) ; jacobLogit
%

% Author: Laura Grigolon and Frank Verboven
% August 2012;

%% Contraction Mapping
[d,iterat]   = DeltaCalculation(theta,Data);

%% Contraction Mapping
% if max(isnan(d)) == 1 || iterat>2500 || isreal(d)==0
if max(isnan(d)) == 1 || isreal(d)==0
 f  = 1e10;
else
%% GMM
if max(isnan(d)) == 1
 f  = 1e10;	  
else 
% Use relationship between linear and non-linear parameters from step 4,
% resulting from the FOC's
    beta  = Data.invxzwzx*(Data.xzwz*d);
    if beta(1)>=-.01
        beta(1)=-.01;
    end
%     if beta(1)>=-.0
%         beta(1)=-.0;
%     end
    % error term
    csi     = d - Data.Xexo*beta;
	f       = csi'*Data.Z*Data.W*Data.Z'*csi;
    save beta beta;
end
%% Gradient
% if max(isnan(d)) == 1 || iterat>2500 || isreal(d)==0    
if max(isnan(d)) == 1 ||  isreal(d)==0    
    % isreal(d)==0 can happen because of the negative weights - this is a problem of quadrature methods
    df=1e10;
else
        temp    = jacob(theta,d,Data);
        df      = 2*temp'*Data.Z*Data.W*Data.Z'*csi;
end

end

