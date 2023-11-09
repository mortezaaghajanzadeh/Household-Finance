clear all
clc
close all

data_lc = csvread('lc_profile.csv',1,0);

% Parameters
Tw = 40;
Tr = 20;
T = Tw+Tr;

beta = 0.94;
gamma = 1.5;
delta = 0;

rw = 0.04;           % Interest rate during working age
rr = 0.04;           % Interest rate during retirement

% Standard deviation of shocks
sigma_e = sqrt(0.015);

% Simulation variables
N = 5000;
mu_a0 = -2.5; % Distribution of initial assets
sigma_a0 = 2;

% Stochastic matricies
A0 = max( exp( mu_a0 + sigma_a0 * min(randn(N,1),3) ) , 0 );  % Initial asset choice
N_s = sigma_e*randn(N,Tw);
unemp = rand(N,Tw);
death = rand(N,T);

%rho_iter = 20; % Set to one for baseline case
%RHO_STORE = linspace(0.7,0.995,rho_iter);

rho_iter = 10; % Set to one for baseline case
RHO_STORE = linspace(0.7,0.95,rho_iter);

PHI_STORE = nan(1,rho_iter);

% Loop over different persistence parameters
% for i0=1:rho_iter
rho_iter = 1;
i0 =1;
if rho_iter==1
rho = 0.97;  % Persistence of income in baseline
else
rho = RHO_STORE(i0);
end
disp(['iter: ', num2str(i0), ' Rho: ',num2str(rho)])

% Income
ny = 5;      % size of income grid
floden_w = 0.5 + rho/4;
sigmaY = sigma_e/sqrt(1-rho^2);
floden_sigma = floden_w*sigma_e + (1-floden_w)*sigmaY;
                       
[ly,Py] = tauchenHussey(ny,0,rho,sigma_e,floden_sigma);  % AR1 process for income
y = exp(ly);

% Transitory component
nu = 2;
pu = 0.01;
Pu = [pu,1-pu];
u = [0.4,1];

P = kron(Py,Pu); % Combine income transition matrix

% Life cycle parameters
%p_dth = ones(1,T); % No death
p_dth = [ones(1,Tw),linspace(1,0.92,Tr)];
%g = cumprod( [1, 1.01 * ones(1,Tw-1)] );
g = data_lc(:,2);

pens = 0.6 * g(Tw);
vmin = -1.e10;

% asset grid
brrw = 0;
na = 10;
amin = 0;
amax = 80;
A = linspace(amin,amax,na)';

Vr = zeros(na+1,Tr);
Cr = zeros(na+1,Tr);
Mr = zeros(na+1,Tr);

% Set up period T from A=0
Mr(2:end,Tr) = A+0.01;
Cr(2:end,Tr) = Mr(2:end,Tr);
Vr(2:end,Tr) = (Cr(2:end,Tr).^(1-gamma)-1)./(1-gamma);
Vr(1,:) = vmin;
% retirment period
for t=Tr-1:-1:1
Mp = A *(1+rr) + pens - rr*brrw;  % Cash on hand tomorrow
Cp = interp1(Mr(:,t+1),Cr(:,t+1),Mp,'linear','extrap');  % interpolate consumption

% Construct tomorrows value and tomorrows derivative wrt assets
EV =  beta * p_dth(t+Tw+1) * interp1(Mr(:,t+1),Vr(:,t+1),Mp,'linear','extrap');
dV = beta * p_dth(t+Tw+1) * Cp.^(-gamma)*(1+rr);

Cr(2:end,t) = dV.^(-1/gamma);  % Use FOC to find consumption policy
Mr(2:end,t) = Cr(2:end,t)+A;   % Implied cash on hand
Vr(2:end,t) = (Cr(2:end,t).^(1-gamma)-1)./(1-gamma) + EV;
end

% Working period
Vw = zeros(ny*(na+1),Tw);
Cw = zeros(ny*(na+1),Tw);
Mw = zeros(ny*(na+1),Tw);

for i=1:ny
Vw((i-1)*na+1*i,:) = vmin;
end
for t=Tw:-1:1
if t==Tw
% Last period working (essentially retirement problem, could be written
% more efficiently by combining this into loop above
Mp = A*(1+rr) + pens - rr*brrw;
Cp = interp1(Mr(:,1),Cr(:,1),Mp,'linear','extrap');
% Construct tomorrows value and tomorrows derivative wrt assets (no
% expectation)
EV = beta * p_dth(t+1) * interp1(Mr(:,1),Vr(:,1),Mp,'linear','extrap');
dV = beta * p_dth(t+1) * Cp.^(-gamma)*(1+rr);

    for i=1:ny
    index = (i-1)*na + 1 + i*1: i*na + i*1;
    Cw(index,t) = dV.^(-1/gamma);
    Mw(index,t) = Cw(index,t)+A;
    Vw(index,t) = (Cw(index,t).^(1-gamma)-1)./(1-gamma) - delta + EV;
    end 

else
Cp = zeros(na,ny*nu);
Vp = zeros(na,ny*nu);
%dVp = zeros(na,ny);
    for i=1:ny
        for j=1:nu  % Loop over transitory shocks
        Mp = A *(1+rw) + y(i)*u(j)*g(t+1) - rw*brrw;  % Implied cash on hand tomorrow
        index = (i-1)*na + i*1: i*na + i*1; % (i-1)*na + 1 + i*1: i*na + i*1;
        index2 = (i-1)*nu + j;
        Cp(:,index2) = interp1(Mw(index,t+1),Cw(index,t+1),Mp,'linear','extrap');
        Vp(:,index2) = interp1(Mw(index,t+1),Vw(index,t+1),Mp,'linear','extrap');    
        end
    end
    
dVp = Cp.^(-gamma);  

% Construct tomorrows value and tomorrows derivative wrt assets (no
% expectation)
EV = beta * p_dth(t+1) * Vp * P';
dV = beta * p_dth(t+1) * dVp * P' * (1+rw);  % RHS of Euler equation
    for i=1:ny
    index = (i-1)*na + 1 + i*1: i*na + i*1;
    Cw(index,t) = dV(:,i).^(-1/gamma);  % Use FOC to find consumption
    Mw(index,t) = Cw(index,t)+A;  % Implied cash on hand
    Vw(index,t) = (Cw(index,t).^(-gamma)-1)/(1-gamma) - delta + EV(:,i);
    end

end

end
disp('Found value function')
%%
% Simulation
A_s = zeros(N,T);
Y_s = zeros(N,T);
U_s = zeros(N,T);
income_s = zeros(N,T);
Yi_s = zeros(N,T);
M_s = zeros(N,T);
C_s = zeros(N,T);
AL_s = ones(N,T);

% Initialise values
A_s(:,1) = A0;
Y_s(:,1) = 1;
U_s(:,1) = u(2);
Yi_s(:,1) = dsearchn(y,Y_s(:,1));

for t=1:T
for i=1:N

if t<=Tw
income_s(i,t) = Y_s(i,t)*U_s(i,t)*g(t);
M_s(i,t) = A_s(i,t)*(1+rw) + income_s(i,t) - rw*brrw; % cash in hand today
% M_s(i,t) = max(min( M_s(i,t), max(Mw(:,t))) ,0);
index = (Yi_s(i,t)-1)*na + Yi_s(i,t)*1: ...
    Yi_s(i,t)*na + Yi_s(i,t)*1;  % index to pick out set of asset choices | on (nearest) income state

C_s(i,t) = interp1(Mw(index,t),Cw(index,t),M_s(i,t),'linear','extrap');  % Consumption choice given state
    % Next period labor income
    if t<Tw        
        if unemp(i,t+1)<Pu(1)
        U_s(i,t+1) = u(1);
        else
        U_s(i,t+1) = u(2);
        end
        Y_s(i,t+1) = exp( rho*log(Y_s(i,t))+N_s(i,t+1) );
        Yi_s(i,t+1) =  dsearchn(y,Y_s(i,t+1));
    end
else
    % Retired households
income_s(i,t) = pens;
M_s(i,t) = A_s(i,t)*(1+rr) + income_s(i,t) - rw*brrw;
% M_s(i,t) = max(min( M_s(i,t), max(Mr(:,t-Tw))) ,0);
C_s(i,t) = interp1(Mr(:,t-Tw),Cr(:,t-Tw),M_s(i,t),'linear','extrap');
end

if t<T
% update assets next period
A_s(i,t+1) =  M_s(i,t)-C_s(i,t);
end

    % Next period's alive or dead
    if t<T
    AL_s(i,t+1) = AL_s(i,t);
    if death(i,t+1)>p_dth(t+1)
    AL_s(i,t+1) = 0;
    end

    end

end

if t==Tw
disp('Completed simulation for working age')
end
end
%%
ALagg = sum(AL_s);  % Population | age
maxA_s = max(A_s);
disp(['Max A choice: ',num2str(max(maxA_s))])

logC_s = log(max(C_s,eps));

% Find aggregates lifecycle paths
Clc = sum(C_s.*AL_s)./ALagg ;
Mlc = sum(M_s.*AL_s)./ALagg ;
Alc= sum(A_s-brrw.*AL_s)./ALagg ;
incomelc = sum(income_s.*AL_s)./ALagg ;

% Aggregate moments 
Agg_Cw = sum(reshape(C_s(:,1:Tw),N*Tw,1).*reshape(AL_s(:,1:Tw),N*Tw,1))/(sum(reshape(AL_s(:,1:Tw),N*Tw,1)));
Agg_Cr = sum(reshape(C_s(:,Tw+1:end),N*Tr,1).*reshape(AL_s(:,Tw+1:end),N*Tr,1))/(sum(reshape(AL_s(:,Tw+1:end),N*Tr,1)));
Agg_C = sum(C_s(:).*AL_s(:))/(sum(AL_s(:)));
Agg_sav = sum((A_s(:)-brrw).*AL_s(:))/(sum(AL_s(:)));

disp(['Av consumption (working): ',num2str(Agg_Cw)])
disp(['Av consumption (retired): ',num2str(Agg_Cr)])
disp(['Agg consumption: ',num2str(Agg_C)])
disp(['Agg savings: ',num2str(Agg_sav)])

% Log.  Consumption growth
logC_s = log(max(C_s,eps));
logClc = sum(logC_s.*AL_s)./ALagg ;
max_logClc = max(logClc);
grwth_logC = max_logClc - logClc(1);
disp(['Growth of log consumption: ',num2str(grwth_logC) ])

logc_s = logC_s - repmat(logClc,N,1);

% Log Income growth
logincome_s = log(max(income_s,eps));
logincomelc = sum(logincome_s.*AL_s)./ALagg ;

max_logincomelc = max(logincomelc);
grwth_logincome = max_logincomelc - logincomelc(1);
disp(['Growth of log income: ',num2str(grwth_logincome)])

% Cohort variance
var_coh_c = var(logC_s(:,1:Tw));
var_coh_inc = var(logincome_s(:,1:Tw));

% BPP insurance
dlogC_s = logC_s(:,2:end) - logC_s(:,1:end-1);
dlogc_s = logc_s(:,2:end) - logc_s(:,1:end-1);

c_temp=reshape( dlogC_s(:,1:Tw-1),N*(Tw-1),1);
n_temp=reshape( N_s(:,2:Tw),N*(Tw-1),1);
cov_cn = cov(c_temp,n_temp);
phi = 1-cov_cn(1,2)/cov_cn(2,2);
PHI_STORE(i0) = phi;

disp(' ')
disp(['var of consumption growth: ',num2str(cov_cn(1,1))]);
disp(['var of perm shocks: ',num2str(cov_cn(2,2))]);
disp(['cov of consumption growth,shocks: ',num2str(cov_cn(1,2))]);
disp(['Partial insurance coefficient: ', num2str(phi)]);
disp(' ')

% end


%%
close all
% Plot Value function
figure(1)
plot( Mr(1:round(na/2),:),Vr(1:round(na/2),:) );
% plot( Mr(1:round(na),:),Vr(1:round(na),:) );
ylabel('V')
%%
figure(2)
plot( Mw(1:round(na/2),:), Vw(1:round(na/2),:) );
ylabel('V')
%%

age_slt_r = round(linspace(1,Tr,5));
age_slt_w = round(linspace(1,Tw,10));
figure(3)
hold on
plot( Mr(1:round(na/2),age_slt_r),Vr(1:round(na/2),age_slt_r),'--' );
plot( Mw(1:round(na/2),age_slt_w), Vw(1:round(na/2),age_slt_w ) );
hold off
ylabel('V')
%%
close all
% Plot policy functions
figure(4)
plot( Mr(1:round(na/2),:),Cr(1:round(na/2),:) );
text(Mr(round(na/4),1),Cr(round(na/4),1),'rt age=1')
text(Mr(round(na/4),end),Cr(round(na/4),end),'rt age=T^r')
ylabel('C')
%% 
figure(5)
plot( Mw(1:round(na),:), Cw(1:round(na),:) );
text(Mw(round(na/4),1),Cw(round(na/4),1),'age=1')
text(Mw(round(na/4),end),Cw(round(na/4),end),'rt age=T^w')
ylabel('C')
%%
figure(6)
hold on
plot( Mr(1:round(na/2),age_slt_r), Cr(1:round(na/2),age_slt_r),'--' );
plot( Mw(1:round(na/2),age_slt_w), Cw(1:round(na/2),age_slt_w ) );
text(Mr(round(na/4),1),Cr(round(na/4),1),'rt age=1')
text(Mr(round(na/4),end),Cr(round(na/4),end),'rt age=T^r')
text(Mw(round(na/4),1),Cw(round(na/4),1),'age=1')
text(Mw(round(na/4),end),Cw(round(na/4),end),'rt age=T^w')
hold off
ylabel('C')

% Plot lifecycle patterns
C_s_plot = C_s;
C_s_plot(AL_s==0)=nan;
income_s_plot = income_s;
income_s_plot(AL_s==0)=nan;

figure(7)
subplot(2,2,1)
plot(1:T,Clc,'LineWidth',3)
title('consumption')

subplot(2,2,2)
plot(1:T,incomelc,'LineWidth',3)
title('income')

subplot(2,2,3)
plot(1:T,Alc,'LineWidth',3)
title('savings')

subplot(2,2,4)
% plot(1:T,Mlc,'LineWidth',3)
% title('cash on hand')
plot(1:T,ALagg/N,'LineWidth',3)
title('% alive')

figure(8)
subplot(2,1,1)
hold on
plot(1:T,C_s_plot,'-m')
plot(1:T,Clc,'LineWidth',3)
hold off
title('consumption')

subplot(2,1,2)
hold on
plot(1:T,income_s_plot,'-m')
plot(1:T,incomelc,'LineWidth',3)
hold off
title('income')

% Plot cohort variance
figure(9)
subplot(1,2,1);
hold on
plot(1:Tw,var_coh_inc,'+-k','LineWidth',2)
hold off
title('var log(y)')
yl = ylim;

subplot(1,2,2)
hold on
plot(1:Tw,var_coh_c,'o-b','LineWidth',2)
hold off
title('var log(c)')
ylim(yl)

% Plot insurance parameter for different income persistence
if rho_iter>1
figure(9)
plot(RHO_STORE,1-PHI_STORE,'LineWidth',3)
xlabel('persistence: \rho')
ylabel('insurance: 1-\phi')
end
