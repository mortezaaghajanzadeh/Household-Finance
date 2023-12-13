%%

clear
close all
clc



useGPU = 0;

%% training set

% nSimul = 100000;
nSimul = 1000;

Allp = haltonset(3,'Skip',1);  Allp = scramble(Allp,'RR2');  Allp = net(Allp,1000) ;
LB = [.8 1.01 0]; UB = [.99 5 .01];
Allp = LB + Allp.*(UB-LB);


%% train
for thiS =  1:size(Allp,1)

disp(thiS)
disp(Allp(thiS,:))

m = simul(Allp(thiS,:),nSimul,useGPU);



disp(m)

Allm(thiS,:) = m;

disp('------------------------')

end




% %%
% 
% clc
% tMoments = true(1,7); %%% just a test, include all
% 
% X = makeitpoly(Allp,2);
% % Y = Allm(1:400,:);
% Y = Allm;
% targetm = readtable('data_moments.txt'); targetm = targetm{:,:};
% fitW = 1./(sum((Y(:,tMoments) - targetm(tMoments)).^2,2)); fitW=fitW/max(fitW);
% betaA = (X'*(fitW.*X))\(X'*(fitW.*Y(:,tMoments)));
% err = Y(:,tMoments) - X*betaA;
% disp(1 - (std(err,fitW,1)./std(Y(:,tMoments),fitW,1)).^2)
% 
% approxloss = @(p) sum((makeitpoly(p,2)*betaA - targetm(tMoments)).^2);
% optimoptions = optimset('MaxFunEvals',1000,'Display','on') ; 
% 
% for ttest = 1:25
%     rng shuffle
%     x0 = LB + rand(1,3).*(UB-LB);
%     [x(ttest,:),fval(ttest,:)] = fminsearchbnd(approxloss,x0,LB,UB,optimoptions);
% %     pause(5)
% end
% [~,idxx] = min(fval);
% x = x(idxx,:)
% 
% truem = simul(x,nSimul,useGPU);
% disp(truem)
% disp(targetm)


