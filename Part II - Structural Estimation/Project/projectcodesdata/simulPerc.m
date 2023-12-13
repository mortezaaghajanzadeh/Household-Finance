%%

clear
close all
clc

useGPU = 0;

%% Parameters

%%% externally calibrated
calibr.mubar = 0.05;
calibr.sigma_mu = 0.15;
calibr.rho = 0;
calibr.rf = .01;
calibr.g = 0.01;
calibr.b = 20;
calibr.T = 45;


%%% foundamentals [to be estimated, random numbers here: as an example]
param.beta = 0.92;
param.gamma = 5;
param.phi = .005;





%% Options

% ---
R.N = 51;

w.N= 300;
w.min = -3;
w.max = 3;

% --- 
c.N = 199;
c.min = .01;
c.max = .99;

pi.N = 101;
pi.min = 0;
pi.max = 1;


%% make grids, exogeneous/state

sigma_shock = calibr.sigma_mu * (1-calibr.rho^2)^0.5;
[R.transR,R.mu] = tauchen(calibr.rho , 0 , sigma_shock^2 , R.N , 2);
clear sigma_shock


%%% just a check that tauchen works nice...
[Vp,Dp] = eig(R.transR');
INDXp = abs(diag(Dp) - 1) < 1e-3;
PiR = Vp(:,INDXp);
PiR = PiR'/sum(PiR);
disp("mean rdev, true/simul: " + 0 + " / " + PiR*R.mu)
disp("std r, true/simul: " + calibr.sigma_mu + " / " + sqrt(PiR*(R.mu.^2)))
disp("rho r, true/simul: " + calibr.rho + " / " + (PiR.*R.mu')*R.transR*R.mu/(PiR*(R.mu.^2)))
clear Vp Dp INDXp Vp PiR
%%%


R.Grid = reshape(exp(R.mu + calibr.mubar - calibr.sigma_mu^2/2 + calibr.rf) , [1,1,R.N]);
calibr.Rf = exp(calibr.rf);

w.Grid = reshape(linspace(w.min,w.max,w.N) , [w.N,1]);




%% make grid, choice vars
c.Grid =  linspace(c.min,c.max,c.N);
pi.Grid = linspace(pi.min,pi.max,pi.N);

c.Grid = repmat(reshape(c.Grid,[1 c.N]) , [pi.N 1]); c.Grid = c.Grid(:); c.Grid = reshape(c.Grid,[1 c.N*pi.N 1]);
pi.Grid = repmat(reshape(pi.Grid,[pi.N 1]) , [1 c.N]); pi.Grid = pi.Grid(:); pi.Grid = reshape(pi.Grid,[1 c.N*pi.N 1]);



%% make grid, value function
V = nan(w.N,R.N,calibr.T+1);
optC = nan(w.N,R.N,calibr.T);
optPi = nan(w.N,R.N,calibr.T);


%% transfer to GPU

switch useGPU

    case true
    
    V = gpuArray(V);
    R.Grid = gpuArray(R.Grid); R.transR = gpuArray(R.transR);
    c.Grid = gpuArray(c.Grid);
    pi.Grid = gpuArray(pi.Grid);
    w.Grid = gpuArray(w.Grid);

    case false

        warning('computation is not delivered to GPU... takes much longer time!')

end






%% backward iteration on the value function


%%% retirement value wrt wealth: 
V(:,:,end) = repmat(calibr.b*exp(w.Grid).^(1-param.gamma)./(1-param.gamma) , [1 R.N]);




%%% backward over the working age
tic
for t=calibr.T:-1:1
clc
disp("age in the backward optimization: " + t)

Wp = ( exp(w.Grid) + exp(calibr.g*t) - param.phi.*(pi.Grid>0) ).*(1 - c.Grid).*(pi.Grid.*R.Grid + (1-pi.Grid)*calibr.Rf); 
Wp(Wp<=0) = nan;

Rp = repmat(R.Grid,[w.N,c.N*pi.N,1]);
     
% tic
    EV=reshape(reshape(param.beta*-exp(interpn(w.Grid(:),R.Grid(:),log(-V(:,:,t+1)),log(Wp),Rp)),[],R.N)*R.transR',[w.N , c.N*pi.N , R.N]);     
% toc        
    util = (c.Grid.*(exp(w.Grid) + exp(calibr.g*t) - param.phi.*(pi.Grid>0))).^(1-param.gamma)./(1-param.gamma);

    [Vt , idxx] =  max(util + EV , [] , 2 );
    V(:,:,t) = squeeze(Vt);
    idxx = squeeze(idxx);

    optC(:,:,t) = c.Grid(idxx);
    optPi(:,:,t) = pi.Grid(idxx);


end
toc



%%

h = figure('Name','policy function');
subplot(1,2,1)
plot(squeeze(optPi(find(w.Grid>0,1),(R.N+1)/2,:)))
subplot(1,2,2)
plot(squeeze(optC(find(w.Grid>0,1),(R.N+1)/2,:)))
h.Position(3) = 2*h.Position(3);


