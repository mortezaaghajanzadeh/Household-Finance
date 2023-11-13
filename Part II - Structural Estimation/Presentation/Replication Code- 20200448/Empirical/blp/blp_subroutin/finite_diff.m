e=10^(-3);

[~,~,sj,~,~, ~,~,~,~] = ShareCalculation(thetaRCNL,deltRCNL,Data);
sj1=sj;
[~,~,sj,~,~, ~,~,~,~]  = ShareCalculation(thetaRCNL+[0;e],deltRCNL,Data);
derShareRC_approxi=(sj-sj1)/e;


[~,~,sj,~,~, ~,~,~,~] = ShareCalculation(thetaRCNL,deltRCNL,Data);
sj1=sj;
[~,~,sj,~,~, ~,~,~,~]  = ShareCalculation(thetaRCNL+[e;0],deltRCNL,Data);
derShareSigseg_approxi=(sj-sj1)/e;


[~,~,sj,~,~, ~,~,~,~] = ShareCalculation(thetaRCNL,deltRCNL,Data);
sj1=sj;
[~,~,sj,~,~, ~,~,~,~]  = ShareCalculation(thetaRCNL,deltRCNL+[e;zeros(nobs-1,1)],Data);
derShareDelta_approxi=(sj-sj1)/e;


derShareDelta_approxigmm_bar(theta20RCNL,Data);