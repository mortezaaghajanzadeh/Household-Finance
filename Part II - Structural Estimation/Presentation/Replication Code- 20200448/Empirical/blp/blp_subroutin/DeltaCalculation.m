function [delta_hat,i] = DeltaCalculation(theta,Data)

% DeltaNL - Contraction Mapping for Random Coefficient Nested Logit
% Works like Logit apart for the weighting of the difference between
% calculated and observed shares
% Syntax:  f = deltaNL(thetaNL,Data)
%
% Inputs:
%    input1 - Non linear parameters
%    input2 - Data
%
% Outputs:
%    Mean Value Delta
%
% MAT-files required: mvaloldNL
% Subfunctions: NLShareCalculation
%

% Author: Laura Grigolon and Frank Verboven
% August 2012;

load mvalold
k   = 100;
km  = 1e-14;
i   = 0;
% if max(isnan(delta0))  == 1 || isreal(delta0) == 0
if max(isnan(delta0))  == 1
    deltastart         =  zeros(Data.nobs,1);
else
    deltastart         = delta0;
end

% Unpack
logobsshare              = Data.logobsshare;

while k > km
    % Market Share
    [~,sh] = ShareCalculation(theta,deltastart,Data);
    delta_hat  = deltastart+(logobsshare-log(sh));
    if max(isnan(delta_hat)) == 1
        disp('No Convergence - delta calculation failed: OVERFLOW')
        break
    end
    
    if isreal(delta_hat) == 0
        break
    end
    
    i = i + 1;
    if i > 3500
        %         disp('No Convergence - delta convergence failed')
        break
    end
        k           = max(abs(delta_hat-deltastart));
        deltastart  = delta_hat;
   
end

disp(['sigma1: ' num2str(theta(1)) ' sigma2: ' num2str(theta(end)) ...
    ' # of iterations for delta convergence:  ' num2str(i) '; max(abs(delta1-deltastart)):  ' num2str(k)])

if isreal(delta_hat)   == 0
    delta0          =  zeros(Data.nobs,1);
else
    delta0 = delta_hat;
end
save mvalold delta0

