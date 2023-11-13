function [se,varcovar]=seblp(bet_hat,theta_hat,Data)

% seblp
% Computes standard errors for all parameters
% Written by Mathias Reynaert (2013)

Z = Data.Z;
W = Data.W;
Xexo = Data.Xexo;
% theta2=parameters(nlin+1:nlin+nrc,1);
deltaopt = DeltaCalculation(theta_hat,Data);
derdel=jacob(theta_hat,deltaopt,Data);
derksi=[-Xexo derdel]'*Z;
vv=inv(derksi*W*derksi');
xi=deltaopt-Xexo*bet_hat;

% covg = zeros(size(Z,2));
% for ii =1:length(Z),
%     covg = covg + Z(ii,:)'*Z(ii,:)*(xi(ii)^2);
% end

covg=zeros(size(Z,2));
covg=Z'*(bsxfun(@times,Z,xi.^2));

varcovar = vv*derksi*W*covg*W*derksi'*vv;
se=sqrt(diag(varcovar));
end 