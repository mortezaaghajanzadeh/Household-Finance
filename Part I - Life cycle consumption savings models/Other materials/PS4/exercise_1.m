clear variables
addpath functions

%% 0. MODEL PARAMETERS
beta = 0.94;  % Consumer's discount factor
gamma = 1.5;  % Consumer's coefficient of RRA
psi = 0.6;    % Pension replacement rate
r = 0.04;     % Interest rate
Nw = 40;      % Number of working periods
Nr = 20;      % Number of retired periods

% Survival probabilities
ps = [repelem(1,1,Nw),linspace(1,0.92,Nr),0];  % Linear decrease

% Discretization parameters
amin = 0;    % Borrowing constraint
amax = 80;   % Upper bound on asset holdings
nA = 1000;   % Asset grid size
nYP = 5;     % Permanent income grid size
nYT = 2;     % Transitory income grid size

% Permanent income process parameters, AR(1) model
rho = 0.97;      % AR(1) coef on permanent income
sigeps = 0.015;  % std. dev. of innovations

% Transitory income process parameters
yt = [0.4;1];      % States
Pt = [0.01;0.99];  % Probabilities

% Simulation parameters
I = 5000;  % Number of simulated households
a0_m = -2.5;
a0_s = 2;
a0max = exp(a0_m + 3*a0_s);

% Other parameters
baseAge = 25;
G = readmatrix("data/lc_profile.csv");
G = G(:,2);
method = 'linear';  % Interpolation method for EGM

%% 1. SET-UP
N = Nw + Nr;

% Discretize log-income process using Tauchen-Hussey (1991)
w = 0.5 + rho/4;  % Floden (2008) weights
baseSigma = w * sigeps + (1-w) * sigeps/sqrt(1-rho^2);
[s,Pp] = tauchenHussey(nYP,0,rho,sigeps,baseSigma);
yp = exp(s);  % nY-by-1 vector
clear loga w baseSigma s

% Construct transition matrix over both types of income shocks
P = repmat(kron(Pt',Pp),nYT,1);

% Define grids
a = linspace(amin,amax,nA)';  % Uniform grid
Xw = cell(1,Nw);  % Cash-on-hand grid during working life
y = reshape(repmat(yt',nYP,1),nYP*nYT,1) .* repmat(yp,nYT,1);
for j = 1:Nw
    Xw{j} = (1+r) * a + G(j) * y';
end
Xr = (1+r) * a + psi*G(Nw);  % Cash-on-hand during retirement

%% 2. SOLVE HOUSEHOLD PROBLEM
pars = {a,beta,gamma,r,N,ps};
[V,c] = solveLifeCycleHH(Xw,Xr,P,pars,method);

%% 3. DIAGNOSTICS, HOUSEHOLD PROBLEM
idy = 8;  % "Normal" income state
idj = 1:10:N;
plX = zeros(nA,size(idj,2));
plY = zeros(nA,size(idj,2));
for k = 1:size(idj,2)
    if idj(k) > Nw
        plX(:,k) = Xr;
        plY(:,k) = c{idj(k)} / (psi*G(Nw));
    else
        plX(:,k) = Xw{idj(k)}(:,idy);
        plY(:,k) = c{idj(k)}(:,idy) / (G(idj(k)) * y(idy));
    end
end

fig1 = figure();
plot(plX,plY,LineWidth=1)

clear idy idj plX plY

%% 4. SIMULATION
a0 = min(exp(a0_m + a0_s*randn(I,1)),a0max);
sT = [zeros(I,1),rand(I,Nw-1) < Pt(1)];  % 1 if transitory income shock
logPi = zeros(I,Nw);     % log(P) over time ~ AR(1) process
ai = [a0,zeros(I,N-1)];  % Assets
xi = zeros(I,N);         % Cash-on-hand
ci = zeros(I,N);         % Consumption
di = zeros(I,N);         % 0 if alive, 1 if dead

% Simulate permanent income process
eps = sigeps * randn(I,N-1);
for j = 2:Nw
    logPi(:,j) = rho * logPi(:,j-1) + eps(:,j-1);
end

% Simulate working life (2 -> Nw)
for j = 1:Nw
    ui = yt(2) - sT(:,j) * (yt(2) - yt(1));
    yi = exp(logPi(:,j)) .* ui * G(j);
    xi(:,j) = (1+r) * ai(:,j) + yi;
    [~,idY] = min(abs(yi - y'),[],2);  % Closest discrete income state
    for i = 1:I  % Interpolate consumption policy
        ci(i,j) = interp1(Xw{j}(:,idY(i)),c{j}(:,idY(i)),xi(i,j),method,"extrap");
    end
    ci(:,j) = max(min(ci(:,j),xi(:,j)-amin),0);  % Enforce constraints
    ai(:,j+1) = xi(:,j) - ci(:,j);
end
clear ui yi

% Simulate retirement (Nw+1 -> N)
for j = Nw+1:N-1
    di(:,j) = di(:,j-1) | (rand(I,1) > ps(j));  % People die
    xi(:,j) = (1+r) * ai(:,j) + psi*G(Nw);
    ci(:,j) = interp1(Xr,c{j},xi(:,j),method,"extrap");
    ci(:,j) = max(min(ci(:,j),xi(:,j)-amin),0);  % Enforce constraints
    ai(:,j+1) = xi(:,j) - ci(:,j);
end
di(:,N) = di(:,N-1) | (rand(I,1) > ps(N));
xi(:,N) = (1+r) * ai(:,N) + psi*G(Nw);
ci(:,N) = xi(:,N);

yi = xi - (1+r) * ai;  % Compute income
logci = log(ci);
logyi = log(yi);
dlogci = logci(:,2:N) - logci(:,1:N-1);
dlogyi = logyi(:,2:N) - logyi(:,1:N-1);

%% 5. LIFE CYCLE PLOTS (Questions in 1.3)
age = baseAge + (0:N-1);

fig1a = figure();  % Assuming all households stay alive until j=N
yyaxis left
plot(age,mean(ci,1),age,mean(yi,1),LineWidth=1)
xline(baseAge+Nw,':')
yyaxis right
plot(age,mean(ai,1),LineWidth=1)
xlabel('Age')
legend('consumption','income','assets (RHS)')

fig1b = figure();
plot(age(Nw+1:N),mean(1-di(:,Nw+1:N)),LineWidth=1)
xlabel('Age')

fig2 = figure();  % Assuming all households stay alive until j=N
plot(age(2:Nw),mean(dlogci(:,1:Nw-1)), ...
     age(2:Nw),mean(dlogyi(:,1:Nw-1)), ...
     LineWidth=1 ...
     )
ylim([-0.1,0.1])
yline(0,'--')
xlabel('Age')
legend('\Delta log-consumption','\Delta log-income')

fig3 = figure();  % NB: Accounting for deaths makes virtually no difference
plot(age(1:Nw),var(logci(:,1:Nw)), ...
     age(1:Nw),var(logyi(:,1:Nw)), ...
     LineWidth=1 ...
     )
xlabel('Age')
legend('log-consumption','log-income')

saveas(fig1a,'plots/avg_cia.png')
saveas(fig1b,'plots/avg_survival.png')
saveas(fig2,'plots/growth_log_ci.png')
saveas(fig3,'plots/variance_log_ci.png')

%% 6. INSURANCE STUFF
covCe = cov([dlogci(:),eps(:)]);
phi = 1 - covCe(1,2)/sigeps^2;
