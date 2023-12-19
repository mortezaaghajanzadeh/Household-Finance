%%
clear
close all
clc
useGPU = 0;
%%
nSimul = 100000;
%%% Estimated Values
load("estimated_point")
param.beta = x(1);
param.gamma = x(2);
param.phi = x(3);
% m = simul([param.beta param.gamma param.phi],nSimul,useGPU);
%% Load training data
load("Allp.mat")
load('Allm.mat')
training.Allp = Allp;
training.Allm = Allm;

%% Estimation boundry
constrains.LB = [.8 1.01 0];
constrains.UB = [.99 5 .01];

%% Load Moments

tMoments = false(1,7); %%% just a test, include all

for v = [1 3 4 6] % 1-participation rate 3-conditional mean risky share 4-mean wealth 6-mean wealth at retirement
    tMoments(v) = true;
end
moment.tMoments = tMoments;


%% Generate the Nosie
n_noise = 10000;
targetm = readtable('data_moments.txt'); targetm = targetm{:,:};
moment.targetm = targetm;
mean_wealth_noisy_moments = targetm(4) + 0.1 * randn(n_noise,1);


%% Estimations
estimated_values = zeros(n_noise,3);
for i = 1:n_noise
    clc
    disp(i)
    moment.targetm(4) = mean_wealth_noisy_moments(i);
    estimated_values(i,:) = point_estimate(training,moment,constrains);
end
moment.targetm = targetm;
mean(estimated_values), std(estimated_values),var(estimated_values)

%%
[sorted_moments, sort_idx] = sort(mean_wealth_noisy_moments);
sorted_params = estimated_values(sort_idx, :);




labels = ["$\hat{\beta}$" "$\hat{\gamma}$" "$\hat{\phi}$"];
titles = ["$\beta$" "$\gamma$" "$\phi$"];
names = ["beta" "gamma" "phi"];

for i = 1:3
    figure(i)
    plot(sorted_moments, sorted_params(:, i), '-', 'LineWidth', 1.5);
    xlabel('$E[W_{it}]$','Interpreter','latex');
    ylabel(labels(i),'Interpreter','latex','Rotation',1)
    title('Robustness Check for ' + titles(i),'Interpreter','latex');
    grid on;
    fig_name = 'robustness_check_parameter' + names(i) + '.png' ;
    saveas(gcf, fig_name);
end
%% Add new moments









%% Estimation for different subsets

estimated_values_subset = zeros(0,3);
mulat = {};

for j = 1:7
    tempt = nchoosek(1:7,j);
    tMoments = false(1,7);
    for i = 1:length(tempt)
        if j <7
            for v = tempt(i,:)
                tMoments(v) = true;
            end
            moment.tMoments = tMoments;
            estimated_values_subset(length(estimated_values_subset) + 1,:) = point_estimate(training,moment,constrains);
            tMoments = false(1,7);
            clc
            disp(length(estimated_values_subset))
        else
            moment.tMoments = true(1,7);
            estimated_values_subset(length(estimated_values_subset) + 1,:) = point_estimate(training,moment,constrains);
        end
    end
    mulat{j} = tempt;
    clear tempt
    
end









