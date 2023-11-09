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

% Define grids
a = linspace(amin,amax,nA)';  % Uniform grid
Xr = (1+r) * a + psi*G(Nw);  % Cash-on-hand during retirement

% Simulation objects
a0 = min(exp(a0_m + a0_s*randn(I,1)),a0max);
sT = [zeros(I,1),rand(I,Nw-1) < Pt(1)];  % 1 if transitory income shock
logPi = zeros(I,Nw);          % log(P) over time ~ AR(1) process
ai = [a0,zeros(I,N-1)];       % Assets
xi = zeros(I,N);              % Cash-on-hand
ci = zeros(I,N);              % Consumption
di = zeros(I,N);              % 0 if alive, 1 if dead
eps = sigeps * randn(I,N-1);  % log-P shocks

%% 2. EVALUATE INSURANCE FOR INCREASING PERSISTENCE OF P = Y/U
pars = {a,beta,gamma,r,N,ps};
simRho = linspace(0.7,0.995,10);
simPhi = zeros(size(simRho));
tic
for m = 1:size(simRho,2)
    % Discretize permanent income process
    w = 0.5 + rho/4;  % Floden (2008) weights
    baseSigma = w * sigeps + (1-w) * sigeps/sqrt(1-rho^2);
    [s,Pp] = tauchenHussey(nYP,0,rho,sigeps,baseSigma);
    yp = exp(s);  % nY-by-1 vector
    P = repmat(kron(Pt',Pp),nYT,1);
    % Construct working-life cash-on-hand grids
    Xw = cell(1,Nw);  % Cash-on-hand grid during working life
    y = reshape(repmat(yt',nYP,1),nYP*nYT,1) .* repmat(yp,nYT,1);
    for j = 1:Nw
        Xw{j} = (1+r) * a + G(j) * y';
    end
    % Solve HH problem
    [~,c] = solveLifeCycleHH(Xw,Xr,P,pars,method);
    % Simulate HH panel
    for j = 2:Nw  % log-permanent income process
        logPi(:,j) = simRho(m) * logPi(:,j-1) + eps(:,j-1);
    end
    for j = 1:Nw  % Working life
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
    for j = Nw+1:N-1  % Retirement
        xi(:,j) = (1+r) * ai(:,j) + psi*G(Nw);
        ci(:,j) = interp1(Xr,c{j},xi(:,j),method,"extrap");
        ci(:,j) = max(min(ci(:,j),xi(:,j)-amin),0);  % Enforce constraints
        ai(:,j+1) = xi(:,j) - ci(:,j);
    end
    xi(:,N) = (1+r) * ai(:,N) + psi*G(Nw);
    ci(:,N) = xi(:,N);
    % Compute Phi
    logci = log(ci);
    dlogci = logci(:,2:N) - logci(:,1:N-1);
    covCe = cov([dlogci(:),eps(:)]);
    simPhi(m) = 1 - covCe(1,2)/sigeps^2;
end
clear i j m w baseSigma s Pp yp P Xw y c
clear ui yi idY logci dlogci covCe
fprintf("Done simulating. Runtime = %.2f seconds",toc)

%% 3. INSURANCE PLOTS
fig1 = figure();
plot(simRho,1-simPhi,'-o',LineWidth=1)
xlabel("\rho"), ylabel("1-\phi")

saveas(fig1,'plots/sim_phi.png')
