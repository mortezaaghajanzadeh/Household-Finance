function x = point_estimate(training,moment,constrains)
    tMoments = moment.tMoments;
    targetm = moment.targetm;
    X = makeitpoly(training.Allp,2);
    Y = training.Allm;
    fitW = 1./(sum((Y(:,tMoments) - targetm(tMoments)).^2,2)); fitW=fitW/max(fitW);
    betaA = (X'*(fitW.*X))\(X'*(fitW.*Y(:,tMoments)));
    err = Y(:,tMoments) - X*betaA;
%     disp(1 - (std(err,fitW,1)./std(Y(:,tMoments),fitW,1)).^2)
    
    approxloss = @(p) sum((makeitpoly(p,2)*betaA - targetm(tMoments)).^2);
    optimoptions = optimset('MaxFunEvals',1000,'Display','on') ; 
    
    for ttest = 1:25
        rng shuffle
        x0 = constrains.LB + rand(1,3).*(constrains.UB-constrains.LB);
        [x(ttest,:),fval(ttest,:)] = fminsearchbnd(approxloss,x0,constrains.LB,constrains.UB,optimoptions);
    %     pause(5)
    end
    [~,idxx] = min(fval);
    x = x(idxx,:);
end

