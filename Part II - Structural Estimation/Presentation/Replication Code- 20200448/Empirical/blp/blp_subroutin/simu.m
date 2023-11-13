clear;
clc;
%  close all
ncb           = 10;
lambda        = 1;    %lambda=1/(1-SigsegTrue) =sigma_zera
% betatrue      = [-1.4; 1];                   % true mean tastes on constant, X1, dummy for segment
% rc_true       = [3.89; 0] ;                             % true random coefficient
% betatrue      = [-2; 1];                   % true mean tastes on constant, X1, dummy for segment
% rc_true       = [3.7; 0] ;                             % true random coefficient
% betatrue      = [-1.0; 1];                   % true mean tastes on constant, X1, dummy for segment
% rc_true       = [2.9; 0] ;                             % true random coefficient
betatrue      = [-1; 1];                   % true mean tastes on constant, X1, dummy for segment
rc_true       = [1.9; 0] ;                             % true random coefficient
nrc           = size(rc_true,1);
nlin          = 1;
NonLinTheta   = [lambda;rc_true];
nseg          = 1;
prods         = nseg*ncb+1;
nmkt          = 9;
nobs          = nmkt*prods;                       % number of observations
bank_char_input = [-1000;-1-log(ncb);0];

[qv qweight]=nwspgr('KPU', 1, 10);
nodes=length(qweight);
qv=qv(:);
mqv=sum(qv.*qweight,1);
qv=bsxfun(@minus,qv,mqv);
qweight = qweight';
qv = qv';
Nnodes=length(qv);

% qv=[0:.1:1];
% qv=qv-mean(qv);
% Nnodes=length(qv);
% qweight = repmat(1/Nnodes,Nnodes,1)';
% qweight = repmat(qweight,nobs,1);
% 
% [qv qweight] = nwspgr('KPN', 1, 7);
% qweight = qweight';
% qv=-exp(qv');
% mqv=sum(qv.*qweight,1);
% qv=bsxfun(@minus,qv,mqv);
% Nnodes=length(qv);

% qv = 0:.1:1;
% qweight = pdf(makedist('Triangular','a',0,'b',1,'c',1),qv);
% qweight = qweight/sum(qweight);
% qv=qv-sum(qv.*qweight);
% Nnodes=length(qv);

figure(2);
hold on;
plot(betatrue(1)+qv.*rc_true(1),qweight,'r','LineWidth',2);
hold off;
% 
id_s          = [1;2*ones(nseg*ncb,1)]; % sector: cash, depo, bond
id_s          = repmat(id_s,nmkt,1);
id_g          = ones(prods,1); % group: each bank with three options
id_g          = repmat(id_g,nmkt,1);
id_t          = kron([1:nmkt]',ones(prods,1));
id_i          = [1;ones(ncb,1)];
id_i          = repmat(id_i,nmkt,1);

id_tg = grp2idx(num2str([id_t id_g]));
dum_tg = sparse(dummyvar(id_tg));
dum_t = sparse(dummyvar(id_t));

id_t_TG = grpstats(id_t,id_tg); % _TG means that collapsed by T and G 
dum_t_TG = sparse(dummyvar(id_t_TG)); % to calculate sum of inclusive value of each group in a market, cannot do at observation level

id_t_T = grpstats(id_t,id_t); % _T means that collapsed by T


X1            = zeros(size(id_s));
X1(id_s==1)   = bank_char_input(1); % cash
X1(id_s==2)   = bank_char_input(2); % dep
X1(id_s==3)   = bank_char_input(3); % dep



%  Data Structure
ffrgrid=[1:nmkt]';
ffr=kron(ffrgrid,ones(prods,1));
Data.ffr            = ffr;
Data.nmkt           = nmkt;
Data.nobs           = nobs;
Data.NonLinTheta    = NonLinTheta;
Data.Nnodes         = Nnodes;
Data.nlin           = nlin;
Data.nrc            = nrc;
Data.X1             = X1;
Data.betatrue       = betatrue;
Data.id_s           = id_s;
Data.qweight          = qweight;
Data.id_tg     = id_tg;
Data.dum_tg     = dum_tg;
Data.id_t     = id_t;
Data.dum_t     = dum_t;
Data.id_t_TG     = id_t_TG;
Data.dum_t_TG     = dum_t_TG;
Data.id_t_T = id_t_T;
Data.qv = qv ;
Data.ncb = ncb;
Data.id_i = id_i;

x=ffr(id_s==2);

%%
fun= @(x)foc(x,Data);
options = optimoptions('fsolve','Display','iter');
[price_hat,fval,exitflag] = fsolve(fun,x,options);

price=ffr;
price(id_s==2)=price_hat;
price(price>ffr)=ffr(price>ffr); % ZLB
Xexo            = [price X1];
Xrandom         = [price id_i];

xvu=zeros(nobs,Nnodes*nrc);
for i = 1:nrc
    xvu(:,(i-1)*Nnodes+1:i*Nnodes)=bsxfun(@times,Xrandom(:,i),qv);
end
delta           = Xexo*betatrue;
% Pack new rate vector to calculate market share
Data.xvu             = xvu;
[sij,sijg,sj,~,~, ~,~,~,~] = ShareCalculation(NonLinTheta,delta,Data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pricegrid1 = grpstats(price(id_i==1&id_s==2),id_t(id_i==1&id_s==2),@(x)mean(x,1));
pricegrid2 = grpstats(price(id_i==0&id_s==2),id_t(id_i==0&id_s==2),@(x)mean(x,1));

sj=sum(sj,2); % sj: num column = group

sjgrid1 = grpstats(sj(id_i==1&id_s==2),id_t(id_i==1&id_s==2),@(x)sum(x,1));
sjgrid2 = grpstats(sj(id_i==0&id_s==2),id_t(id_i==0&id_s==2),@(x)sum(x,1));
h=figure(1);
hold on;
yyaxis left
plot(ffrgrid,pricegrid1,'-')
yyaxis right
plot(ffrgrid,sjgrid1,'d-')
hold off;
legend('cb spread','cb dep');
title(['ncb:' num2str(ncb) ' \beta:' num2str(betatrue(1),2) ' \sigma:' num2str(NonLinTheta(2),2) ...
    ' \lambda:' num2str(NonLinTheta(1),2) ' cash:' num2str(bank_char_input(1),2) ' dep:' num2str(bank_char_input(2),2)])

hold on;

% hold off;