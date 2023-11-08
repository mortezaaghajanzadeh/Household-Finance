clear variables
addpath functions
rng(1234);

%% 0. MODEL PARAMETERS
beta = 0.945;  % Consumer's discount factor (baseline = 0.945)
gamma = 2;     % Consumer's coefficient of RRA
psi = 0.8;     % Pension replacement rate (lambda in PS) (baseline = 0.8)
r = 0.02;      % Interest rate
Nw = 40;       % Number of working periods
Nr = 35;       % Number of retired periods

% Survival probabilities (should only be < 1 after retirement)
% ps = [repelem(1,1,Nw),linspace(1,0.92,Nr),0];  % Linear decrease
ps = [repelem(1,1,Nw+Nr),0];  % No death

% Discretization parameters
amin = 0;    % Borrowing constraint
amax = 1000; % Upper bound on asset holdings
Na = 500;    % Asset grid size
Nyp = 15;    % Permanent income grid size
Nyt = 5;     % Transitory income grid size

% Permanent income process parameters (AR(1), Normal)
rho = 1;          % AR(1) coef on permanent income
sigeta = 0.0713;  % Std. dev. of innovations
sigz0 = 0.589;    % Std. dev. of initial shock
yfloor = 4.8;     % Floor on disposable income

% Transitory income process parameters (Normal)
sigomega = 0.1;  % Std. dev.

% Simulation parameters
I = 1000;  % Number of simulated households
a0_s = 2.129;
a0_m = 1.916 - a0_s^2/2;
a0max = exp(a0_m + 3*a0_s);

% Other parameters
baseAge = 25;
G = readmatrix("data/income_profile.csv");
G = G(:,2);
method = 'linear';  % Interpolation method for EGM
file_suffix = 'base';  % 'base', 'psi0', 'hiBt', 'RBe1'

%% 1. SET-UP
N = Nw + Nr;

% Discretize log-income process using Tauchen-Hussey (1991)
[yp,Pp] = tauchenHussey(Nyp,0,rho,sigeta,sigeta);
yp = exp(yp - sigeta^2/2);  % nY-by-1 vector

% Discretize transitory shocks distribution
[yt,Pt] = ghqnorm(Nyt,-sigomega^2/2,sigomega^2);
yt = exp(yt);

% Define grids
loga = linspace(0,log(amax-amin + 1),Na)';
a = exp(loga) - 1 + amin;  % log-linear grid
clear loga
% a = linspace(amin,amax,Na)';

%% 2. SOLVE HOUSEHOLD PROBLEM
genPars = {a,beta,gamma,ps};
incPars = {G,yfloor,psi,yp,Pp,yt,Pt};
[X,V,c] = solveLifeCycleHH(r,Nw,Nr,genPars,incPars,method);

%% 3. DIAGNOSTICS, HOUSEHOLD PROBLEM
% Plot consumption policy across ages, for different permanent income states
ids = [1,15];  % Income state(s)
idj = 1:10:Nw+1;
plX = cell(size(ids));
plc = cell(size(ids));
plV = cell(size(ids));
maxc = zeros(size(ids));  % To align subplots
minV = zeros(size(ids));
for k = 1:size(ids,2)
    tmpX = zeros(Na,size(idj,2));
    tmpc = zeros(Na,size(idj,2));
    tmpV = zeros(Na,size(idj,2));
    for i = 1:size(idj,2)
        tmpX(:,i) = X{idj(i)}(:,ids(k));
        tmpc(:,i) = c{idj(i)}(:,ids(k));
        tmpV(:,i) = V{idj(i)}(:,ids(k));
    end
    plX{k} = tmpX;
    plc{k} = tmpc;
    plV{k} = tmpV;
    maxc(k) = ceil(max(tmpc,[],'all'));
    minV(k) = floor(min(tmpV,[],'all'));
end
maxc = max(maxc);
minV = mean(minV);

fig1a = figure();
for k = 1:size(ids,2)
    subplot(1,2,k)
    plot(plX{k},plc{k},LineWidth=1); hold on
    colororder(parula(size(idj,2)))
    plot([amin,amax],[amin,amax],'k:')
    ylim([0,maxc]); xlim([amin,amax])
    labels = cell(size(idj));
    for i = 1:size(idj,2)
        labels{i} = sprintf('Age %d',baseAge-1 + idj(i));
    end
    legend([string(labels),""],Location="southeast")
    xlabel('Cash-on-hand')
    if k == 1
        ylabel('Consumption')
    end
    title(sprintf('Income state %d',ids(k)))
    hold off
end

fig1b = figure();
for k = 1:size(ids,2)
    subplot(1,2,k)
    plot(plX{k},plV{k},LineWidth=1);
    colororder(parula(size(idj,2)))
    ylim([max(minV,-5),0]); xlim([amin,amax])  % Ad-hoc ylim
    labels = cell(size(idj));
    for i = 1:size(idj,2)
        labels{i} = sprintf('Age %d',baseAge-1 + idj(i));
    end
    legend(string(labels),Location="southeast")
    xlabel('Cash-on-hand')
    if k == 1
        ylabel('Continuation value')
    end
    title(sprintf('Income state %d',ids(k)))
end
clear ids idj plX plc plV maxc minV tmpX tmpc tmpV

saveas(fig1a,sprintf('plots/%s_hh_c_byage.png',file_suffix))
saveas(fig1b,sprintf('plots/%s_hh_V_byage.png',file_suffix))

%% 4.a SIMULATION SHOCKS
% Initial assets
a0 = min(exp(a0_m + a0_s*randn(I,1)),a0max);

% Transitory income state
st = zeros(I,Nw);  % Takes values in 1:Nyt
Ft = cumsum(Pt);
for j = 1:Nw
    [~,st(:,j)] = max(rand(I,1) <= Ft',[],2);
end
clear Ft

% Initial permanent income
z0 = exp(-sigz0^2/2 + sigz0*randn(I,1));
[~,sP0] = min(abs(z0 - yp'),[],2);

% Permanent income state
sp = [sP0,zeros(I,Nw-1)];  % Takes values in 1:Nyp
Fp = cumsum(Pp,2);
for j = 2:Nw
    for i = 1:I
        [~,sp(i,j)] = find(rand() <= Fp(sp(i,j-1),:),1);
    end
end
clear Fp

% Compute income over time
yi = [G' .* yp(sp) .* yt(st), repmat(psi*G(Nw) * yp(sp(:,Nw)),1,Nr)];

% Death
di = zeros(I,N);         % 0 if alive, 1 if dead
for j = 2:N
    di(:,j) = di(:,j-1) | (rand(I,1) > ps(j));  % People die
end

%% 4.b SIMULATION
ai = [a0,zeros(I,N-1)];  % Assets
xi = zeros(I,N);         % Cash-on-hand
ci = zeros(I,N);         % Consumption

% Simulate life cycle
for j = 1:N-1
    xi(:,j) = (1+r) * ai(:,j) + yi(:,j);
    for i = 1:I  % Interpolate consumption policy
        spij = sp(i,min(j,Nw));  % yp-state at j=Nw prevails through retirement
        ci(i,j) = interp1( ...
            X{j}(:,spij), c{j}(:,spij), xi(i,j), ...
            method,"extrap" ...
        );
    end
    ci(:,j) = max(min(ci(:,j),xi(:,j)-amin),0);  % Enforce constraints
    ai(:,j+1) = xi(:,j) - ci(:,j);
end
% Last period (no bequests)
xi(:,N) = (1+r) * ai(:,N) + yi(:,N);
ci(:,N) = xi(:,N);

%% 5. LIFE CYCLE PLOTS (Questions in 1.3)
age = baseAge + (0:N-1);
cmap = parula(4);

% Average life-cycle profile
fig1 = figure();  % Assuming all households stay alive until j=N

subplot(3,1,1:2)
colororder(cmap([1,3],:))
yyaxis left
plot(age,mean(ci,1),age,mean(yi,1),LineWidth=1)
yyaxis right
plot(age,mean(ai,1),LineWidth=1)
xline(baseAge+Nw,':')
legend('Consumption','Income','Assets (RHS)')
xlim([age(1),age(end)])
set(gca,"FontSize",9)

subplot(3,1,3)
colororder(cmap(3,:))
plot(age,mean(1-ci./xi,1),LineWidth=1)
ylim([0,1]); xlim([age(1),age(end)])
xline(baseAge+Nw,':')
xlabel('Age'); ylabel('Saving rate')
set(gca,"FontSize",9)

% Life-cycle profile for selected individuals
qtl = [0.05,0.95];  % Quantile levels on initial assets and income
qtx = quantile(xi(:,1),qtl);
[~,qti] = min(abs(xi(:,1)-qtx),[],1);

fig2 = figure();
for k = 1:2
    ax = subplot(2,1,k);
    colororder(ax,cmap([1,3],:))
    yyaxis left
    plot(age,ci(qti(k),:),age,yi(qti(k),:),LineWidth=1)
    yyaxis right
    plot(age,ai(qti(k),:),LineWidth=1)
    xline(baseAge+Nw,':')
    legend('Consumption','Income','Assets (RHS)')
    xlim([age(1),age(end)])
    title(sprintf("Individual at %d th percentile of x(25)",qtl(k)*100))
end
xlabel('Age')

saveas(fig1,sprintf('plots/%s_avg_lifecycle.png',file_suffix))
saveas(fig2,sprintf('plots/%s_lifecycle_individuals.png',file_suffix))
