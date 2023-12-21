function [m,sPi,sW] = simul(paramvec,nSimul,useGPU,capital_tax)


%% Input Parameters

param.beta = paramvec(1); % could be anything between [.8,.99]
param.gamma = paramvec(2); % could be anything between [1.01,5]
param.phi = paramvec(3); % could be anything between [0,0.01]



%% externally calibrated
calibr.premium = 0.05;
calibr.sharpe_ = 1/3;
calibr.rf = .01;
calibr.g = 0.01;
calibr.b = 20;
calibr.T = 45;
calibr.tau  = capital_tax; %% Capital Income Tax






%% Options


% ---
R.N = 41;

w.N= 600;
w.min = log(.01);
w.max = log(500);

% --- 
c.N = 101;
c.min = .0001;
c.max = .9999;

pi.N = 101;
pi.min = 0;
pi.max = 1;



%% make grids, exogeneous/state

calibr.sigma_mu = (log((calibr.premium/calibr.sharpe_ / (1 + calibr.rf + calibr.premium ))^2 + 1))^.5;
calibr.mean_mu = log(1 + calibr.rf + calibr.premium ) - calibr.sigma_mu^2/2;

[R.transR,R.mu] = tauchen(0 , 0 , calibr.sigma_mu^2 , R.N , 2.5);
% %%% just a check that tauchen works nice...
% [Vp,Dp] = eig(R.transR');
% INDXp = abs(diag(Dp) - 1) < 1e-3;
% PiR = Vp(:,INDXp);
% PiR = PiR'/sum(PiR);
% disp("mean shock, true/discrete: " + 0 + " / " + PiR*R.mu )
% disp("std shock, true/discrete: " + calibr.sigma_mu + " / " + sqrt(PiR*(R.mu.^2)))
% disp("AR(1) shock, true/discrete: " + 0 + " / " + (PiR.*R.mu')*R.transR*R.mu/(PiR*(R.mu.^2)))
% clear Vp Dp INDXp Vp PiR  
% %%%
R.transR = R.transR(1,:); %%% because AR(1), realization is independent of the previous state


R.Grid = reshape(exp(R.mu + calibr.mean_mu) , [1,1,R.N]);
calibr.Rf = 1+calibr.rf;

w.Grid = reshape(linspace(w.min,w.max,w.N) , [w.N,1]);




%% make grid, choice vars
c.Grid =  linspace(c.min,c.max,c.N);
pi.Grid = linspace(pi.min,pi.max,pi.N);

c.Grid = repmat(reshape(c.Grid,[1 c.N]) , [pi.N 1]); c.Grid = c.Grid(:); c.Grid = reshape(c.Grid,[1 c.N*pi.N 1]);
pi.Grid = repmat(reshape(pi.Grid,[pi.N 1]) , [1 c.N]); pi.Grid = pi.Grid(:); pi.Grid = reshape(pi.Grid,[1 c.N*pi.N 1]);



%% make grid, value function
V = nan(w.N,calibr.T+1);
opt_c = nan(w.N,calibr.T);
opt_pi = nan(w.N,calibr.T);



%% transfer to GPU

switch useGPU

    case true
        V = gpuArray(V);
        R.Grid = gpuArray(R.Grid); 
        c.Grid = gpuArray(c.Grid);
        pi.Grid = gpuArray(pi.Grid);
        w.Grid = gpuArray(w.Grid);

    case false
        warning('computation is not delivered to GPU... takes much longer time!')

end






%% backward iteration on the value function


%%% retirement value wrt wealth: 
V(:,end) = calibr.b*exp(w.Grid).^(1-param.gamma)./(1-param.gamma);




%%% backward over the working age
% tic
for t=calibr.T:-1:1


capital_income = (pi.Grid.*R.Grid* (1- calibr.tau ) + (1-pi.Grid)*calibr.Rf * (1- calibr.tau) + calibr.tau);
Wp = (1 - c.Grid).*( exp(w.Grid) + exp(calibr.g*t)*(1 - param.phi.*(pi.Grid>0) ) ).*( capital_income); 
Wp(Wp<=0) = nan;

     
% tic
    bEV=reshape(reshape(param.beta*-exp(interpn(w.Grid(:),log(-V(:,t+1)),log(Wp))),[],R.N)*R.transR',[w.N , c.N*pi.N]);     
% toc        
    util = (c.Grid.*(exp(w.Grid) + exp(calibr.g*t)*(1 - param.phi.*(pi.Grid>0)) )).^(1-param.gamma)./(1-param.gamma);

    [V(:,t) , idxx] =  max( util + bEV , [] , 2 );

    opt_c(:,t) = reshape(c.Grid(idxx),size(w.Grid));
    opt_pi(:,t) = reshape(pi.Grid(idxx),size(w.Grid));


end
% toc



%% simulate 

rng default


calibr.wimean = -1;
calibr.wistd = 1;


sW = nan(nSimul,calibr.T+1);
sPi = nan(nSimul,calibr.T);
sSR = nan(nSimul,calibr.T);
sC = nan(nSimul,calibr.T);

% Rcommonshock = randn(1,calibr.T);
% sR = exp(Rcommonshock*calibr.sigma_mu + calibr.mean_mu);
sR = readmatrix('realized_R.txt');


wishock = randn(nSimul,1);
sW(:,1) = exp(wishock*calibr.wistd + calibr.wimean - calibr.wistd^2/2);

% disp("mean initial wealth: " + mean(sW(:,1)))
% disp("std initial wealth: " + std(sW(:,1)))
% if mean(log(sW(:,1))<w.min+1 | log(sW(:,1))>w.max-2) > 0 
%     warning('w grid range may be narrow!')
% end



for t = 1:calibr.T

%    outrange = mean(log(sW(:,t))>w.max | log(sW(:,t))<w.min)*100;
%     if outrange>0
%         disp("Nan %: " + outrange)
%         warning('out of grid!')
%     end  

   spi = interpn(w.Grid(:),opt_pi(:,t),log(sW(:,t)),'spline'); spi = max(spi,pi.min); spi = min(spi,pi.max);
   sconsmrate = interpn(w.Grid(:),opt_c(:,t),log(sW(:,t)),'spline'); sconsmrate = max(sconsmrate,c.min); sconsmrate = min(sconsmrate,c.max);
   ssavingrate = ( sW(:,t) + exp(calibr.g*t)*(1 - param.phi.*(spi>0)) ).*(1 - sconsmrate)./sW(:,t);

   sW(:,t+1) = sW(:,t).*ssavingrate.*(spi.*sR(:,t) + (1-spi)*calibr.Rf); 
   sPi(:,t) = spi;
   sSR(:,t) = ssavingrate;
   sC(:,t) = sconsmrate.*( sW(:,t) + exp(calibr.g*t)*(1 - param.phi.*(spi>0)) );

end


%% calculate moments

% figure(1)
% plot(sR,'--')
% title('realized r')
% figure(2)
% plot(sort(sW(:,1)),'r--')
% title('initial wealth sorted')
% figure(3)
% plot(mean(sW(:,:),1))
% title('mean wealth')
% figure(4)
% plot(std(sW(:,:),1,1))
% title('std wealth')
% figure(5)
% plot(mean(sPi(:,:),1))
% title('mean participation')
% figure(6)
% plot(mean(sSR(:,:),1))
% title('mean saving rate')
% figure(7)
% plot(mean(sC(:,:),1))
% title('mean consumption')

m(1) = mean(sPi>0,'all');
m(2) = mean(sPi,'all');
m(3) = mean(sPi(sPi>0),'all');
m(4) = mean(sW,'all');
m(5) = std(sW,1,'all');
m(6) = mean(sW(:,end));
m(7) = std(sW(:,end),1);




end

