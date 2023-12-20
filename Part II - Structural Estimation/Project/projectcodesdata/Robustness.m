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
n_noise = 1000;
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

%% Graphs
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
%     title('Robustness Check for ' + titles(i),'Interpreter','latex');
    yline(x(i),'g')
    mean_value = mean(estimated_values(:,i));
    conf_term_2 = 1.96 * (var(estimated_values(:,i)) / length(estimated_values))^0.5;
    yline(mean_value,'r')
    yline(mean_value + conf_term_2 ,'-.r')
    yline(mean_value - conf_term_2 ,'-.r')
    grid on;
    fig_name = 'robustness_check_parameter' + names(i) + '.png' ;
    saveas(gcf, fig_name);
end



%% Add new moments
close all
moment.targetm = targetm;
moment.tMoments = tMoments;

moment_inclusion = {};
inclusion_moments = [2 5 7];
for i = inclusion_moments
    moment.tMoments(i) = true;
    estimated_value = zeros(n_noise,3);
    noisy_moments = targetm(i) + 0.1 * randn(n_noise,1);
    for j= 1:n_noise
        clc
        disp([i j])
        moment.targetm(i) = noisy_moments(j);
        estimated_value(j,:) = point_estimate(training,moment,constrains);
    end
    moment.tMoments = tMoments;
    moment_inclusion{i} = estimated_value;
    clear estimated_value
end
%% Graph
Y.beta = [param.beta];
std.beta = [0];
Y.gamma = [param.gamma];
std.gamma = [0];
Y.phi = [param.phi];
std.phi = [0];

for i = inclusion_moments
    estimated_value = moment_inclusion{i};
    Y.beta =  [Y.beta mean(estimated_value(:,1))];
    std.beta =  [std.beta var([estimated_value(:,1)])^0.5];
    Y.gamma = [Y.gamma mean(estimated_value(:,2))];
    std.gamma = [std.gamma var([estimated_value(:,2)])^0.5];
    Y.phi = [Y.phi mean(estimated_value(:,3))];
    std.phi = [std.phi var([estimated_value(:,3)])^0.5];
end


%%
X = ["Base line" "E[\pi_{it}]" "\sigma[W_{it}]" "\sigma[W_{iT}]"];
x=categorical(X(2:4));
Y_agg = [Y.beta; Y.gamma ;Y.phi];
std_agg = [std.beta; std.gamma; std.phi];
close all

labels = ["$\hat{\beta}$" "$\hat{\gamma}$" "$\hat{\phi}$"];
titles = ["$\beta$" "$\gamma$" "$\phi$"];
names = ["beta" "gamma" "phi"];

for i=1:3
    figure(i)
    y = Y_agg(i,:);
    conf_term = 1.96 * std_agg(i,:)/ (length(std_agg(i,:)) ^ 0.5);
    errorbar(x,y(2:4),conf_term(2:4),"o")
    mean_value = mean(estimated_values(:,i));
    conf_term_2 = 1.96 * (var(estimated_values(:,i)) / length(estimated_values))^0.5;
    yline(mean_value,'r')
    yline(mean_value + conf_term_2 ,'-.r')
    yline(mean_value - conf_term_2 ,'-.r')
    grid on
    ylabel(labels(i),'Interpreter','latex','Rotation',1)
%     title('Moment Inclusion for ' + titles(i),'Interpreter','latex');
    grid on;
    fig_name = 'moment_inclusion_check_parameter' + names(i) + '.png' ;
    saveas(gcf, fig_name);
end





%%




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









