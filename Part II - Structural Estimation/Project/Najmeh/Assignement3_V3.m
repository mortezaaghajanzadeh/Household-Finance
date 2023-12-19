%%

clear
close all
clc



useGPU = 0;

%% training set

nSimul = 100000;

Allp = haltonset(3,'Skip',1);  Allp = scramble(Allp,'RR2');  Allp = net(Allp,1000) ;
LB = [.8 1.01 0]; UB = [.99 5 .01];
Allp = LB + Allp.*(UB-LB);


%% train
% for thiS =  1:size(Allp,1)
% 
% disp(thiS)
% disp(Allp(thiS,:))
% 
% m = simul(Allp(thiS,:),nSimul,useGPU);
% 
% 
% 
% disp(m)
% 
% Allm(thiS,:) = m;
% 
% disp('------------------------')
% 
% end

load('Allm.mat')



%% Approximate SMM
% % All moments
% 
% clc
% tMoments = true(1,7); 
% 
% X = makeitpoly(Allp,2);
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
% x = x(idxx,:) % x is parameteres 
% 
% truem = simul(x,nSimul,useGPU);
% disp(truem)
% disp(targetm)
%% Approximate SMM

% % Selected moments:m1,m3,m4 
% 
% clc
% tMoments = false(1, 7);tMoments([1, 3, 4]) = 1;
% 
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
% 
% x = x(idxx,:) % x is parameteres  
% truem = simul(x,nSimul,useGPU);
% disp(truem)
% disp(targetm)
% 

%% Approximate SMM
% Selected moments:m1,m3,m4,m6 

clc
tMoments = false(1, 7);tMoments([1, 3, 4, 6]) = 1;


X = makeitpoly(Allp,2);
% Y = Allm(1:400,:);
Y = Allm;
targetm = readtable('data_moments.txt'); targetm = targetm{:,:};
fitW = 1./(sum((Y(:,tMoments) - targetm(tMoments)).^2,2)); fitW=fitW/max(fitW);
betaA = (X'*(fitW.*X))\(X'*(fitW.*Y(:,tMoments)));
err = Y(:,tMoments) - X*betaA;
disp(1 - (std(err,fitW,1)./std(Y(:,tMoments),fitW,1)).^2)

approxloss = @(p) sum((makeitpoly(p,2)*betaA - targetm(tMoments)).^2);
optimoptions = optimset('MaxFunEvals',1000,'Display','on') ; 

for ttest = 1:25
    rng shuffle
    x0 = LB + rand(1,3).*(UB-LB);
    [x(ttest,:),fval(ttest,:)] = fminsearchbnd(approxloss,x0,LB,UB,optimoptions);
%     pause(5)
end
[~,idxx] = min(fval);

x_base = x(idxx, :)  % x_base is the base model parameters

% Display base model results
simul_m_base = simul(x_base, nSimul, useGPU,0);
disp("Base Model Results:")
disp(simul_m_base)
disp(targetm)
%% Robustness 1
clc
% Robustness check with noise in mean wealth moment
std_dev_noise = 0.1;  % Standard deviation of the noise
num_simulations = 100; % Number of robustness check iterations

% Initialize arrays to store parameter values for each iteration
robust_params = zeros(num_simulations, 3);
robust_moments = zeros(num_simulations, 1);

for i = 1:num_simulations
    % Introduce noise to the mean wealth moment
    targetm_noisy = targetm;
    targetm_noisy(4) = targetm_noisy(4) + std_dev_noise * randn;
    
    % Robustness check estimation
    fitW_noisy = 1./(sum((Y(:, tMoments) - targetm_noisy(tMoments)).^2, 2));
    fitW_noisy = fitW_noisy / max(fitW_noisy);
    betaA_noisy = (X' * (fitW_noisy .* X)) \ (X' * (fitW_noisy .* Y(:, tMoments)));
    err_noisy = Y(:,tMoments) - X*betaA_noisy;
    disp(1 - (std(err_noisy,fitW_noisy,1)./std(Y(:,tMoments),fitW_noisy,1)).^2)
   approxloss_noisy = @(p) sum((makeitpoly(p,2)*betaA_noisy - targetm_noisy(tMoments)).^2);
   optimoptions_noisy = optimset('MaxFunEvals',1000,'Display','on') ; 

for ttest = 1:25
    rng shuffle
    x0 = LB + rand(1,3).*(UB-LB);
    [x(ttest,:),fval(ttest,:)] = fminsearchbnd(approxloss_noisy,x0,LB,UB,optimoptions);
%     pause(5)
end
[~,idxx] = min(fval);

x_noisy = x(idxx, :)  % x_noisy is the noisy model parameters

% Store parameter values for analysis
robust_params(i, :) = x_noisy;
robust_moments(i) = targetm_noisy(4);
end

% Sort robust_moments and adjust corresponding rows in robust_params
[sorted_moments, sort_idx] = sort(robust_moments);
sorted_params = robust_params(sort_idx, :);
%%
for i = 1:3
    subplot(3, 1, i);
    plot(sorted_moments, sorted_params(:, i), 'o-', 'LineWidth', 1.5);
    xlabel('Noisy Mean Wealth Moment');
    ylabel(['Parameter ' num2str(i)]);
    title(['Robustness Check for Parameter ' num2str(i)]);
    grid on;
      % Save the figure
    fig_name = ['robustness_check_parameter.png'];
    saveas(gcf, fig_name);
end

%% Robustness 2
clc

moments_to_add = [2, 5, 7];
x_moments_add = cell(length(moments_to_add), 1);

% Iterate over parameters
for param_index = 1:3
    figure; % Create a new figure for each parameter

    hold on;

    % Plot the baseline parameter estimate
    baseline_line = plot([0, length(moments_to_add) + 1], [x_base(param_index), x_base(param_index)], 'k--', 'LineWidth', 2);

    % Iterate over moments to add
    for i = 1:length(moments_to_add)
        moment_index = moments_to_add(i);
        tMoments = false(1, 7);
        tMoments([1, 3, 4, 6, moment_index]) = 1;

        X = makeitpoly(Allp, 2);
        Y = Allm;
        targetm = readtable('data_moments.txt'); targetm = targetm{:,:};

        fitW = 1./(sum((Y(:,tMoments) - targetm(tMoments)).^2,2));
        fitW = fitW/max(fitW);

        betaA = (X'*(fitW.*X))\(X'*(fitW.*Y(:,tMoments)));
        err = Y(:,tMoments) - X*betaA;

        disp(1 - (std(err, fitW, 1)./std(Y(:,tMoments), fitW, 1)).^2)

        approxloss = @(p) sum((makeitpoly(p, 2)*betaA - targetm(tMoments)).^2);
        optimoptions = optimset('MaxFunEvals', 1000, 'Display', 'on');

        for ttest = 1:25
            rng shuffle
            x0 = LB + rand(1, 3).*(UB - LB);
            [x(ttest,:), fval(ttest,:)] = fminsearchbnd(approxloss, x0, LB, UB, optimoptions);
        end

        [~, idxx] = min(fval);
        x_moments_add{i} = x(idxx, :);
        
        disp(['x_moments_add ' num2str(moment_index) ':']);
        disp(x_moments_add{i});
        
        % Plot the parameter in x_moments_add{i}
        scatter(i, x_moments_add{i}(param_index), 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
    end

    hold off;

    title(['Parameter ' num2str(param_index) ' vs Additional Moment']);
    xlabel('Moments Added');
    ylabel(['Parameter ' num2str(param_index) ' Values']);

    legend([baseline_line], {'Baseline Parameter Estimate'}, 'Location', 'Best');
    % Modify x-axis ticks and labels
    xticks([1 2 3]);  % Set the x-axis ticks to be at positions 1, 2, and 3
    xticklabels({'mean risky share', 'std wealth', 'std wealth at retirement'});  % Set the labels for each tick position
    xlabel('Additional Moment');  
    % Save the figure
    fig_name = ['parameter_' num2str(param_index) '_vs_moments_added.png'];
    saveas(gcf, fig_name);

    
end

%%
%How calculate J?  nlinfit?
%ci = nlparci(beta,r,"Jacobian",J)

% % Calculate the standard errors using the SMM weight matrix (fitW)
% W = diag(fitW);
% JWJ_inv = inv(J' * W * J);
% 
% % Calculate the standard errors of the parameter estimates
% se = sqrt(diag(JWJ_inv));
% 
% % Compute the 95% confidence intervals
% alpha = 0.05;
% z_critical = norminv(1 - alpha / 2);
% ci = [x_base - z_critical * se, x_base + z_critical * se];

%% Counterfactual
% Counterfactual 1
m_tax = simul(x_base,nSimul,useGPU,0.3);
disp(simul_m_base(6))
disp(m_tax(6))
% Counterfactual 2
x_nocost = x_base;
x_nocost(3) = 0;
m_nocost = simul(x_nocost,nSimul,useGPU,0);
% Counterfactual 3
x_gamma = x_base;
x_gamma(:, 2) = 2 * x_base(:, 2);
m_gamma = simul(x_gamma,nSimul,useGPU,0);

%% Hetorogeneity 

% Type I
disp("The simulation for type I")
nSimul_1 = nSimul/2;
x_type1 = x_base
x_type1(1) = x_base(:, 1) + 0.01;
x_type1(2) = x_base(:, 2) - 0.5;
x_type1(3) = x_base(:, 3) * 0; 
m_1 = simul(x_type1,nSimul_1,useGPU,0.0); % With no capital tax

% Type II
disp("The simulation for type II")
nSimul_2 = nSimul/2;
x_type2 = x_base
x_type2(1) = x_base(:, 1) - 0.01;
x_type2(2) = x_base(:, 2) + 0.5;
x_type2(3) = x_base(:, 3) * 2; 
m_2 = simul(x_type2,nSimul_2,useGPU,0.0);  % With no capital tax

% Aggregate moments
m_star = (m_1 * nSimul_1 + m_2 * nSimul_2) / (nSimul_1 + nSimul_2);

targetm = readtable('data_moments.txt'); targetm = targetm{:,:};
fitW_star = 1./(sum((Y(:,tMoments) - m_star(tMoments)).^2,2)); fitW_star=fitW_star/max(fitW_star);
betaA = (X'*(fitW_star.*X))\(X'*(fitW_star.*Y(:,tMoments)));
err = Y(:,tMoments) - X*betaA;
disp(1 - (std(err,fitW_star,1)./std(Y(:,tMoments),fitW_star,1)).^2)

approxloss = @(p) sum((makeitpoly(p,2)*betaA - m_star(tMoments)).^2);
optimoptions = optimset('MaxFunEvals',1000,'Display','on') ; 

for ttest = 1:25
    rng shuffle
    x0 = LB + rand(1,3).*(UB-LB);
    [x(ttest,:),fval(ttest,:)] = fminsearchbnd(approxloss,x0,LB,UB,optimoptions);
%     pause(5)
end
[~,idxx] = min(fval);

x_star = x(idxx, :)  % x_base is the base model parameters

