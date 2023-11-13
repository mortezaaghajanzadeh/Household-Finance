function data_simu = simulation(beta_hat,theta_hat,bank_char_hat,qv,qweight)

%%
betatrue      = [beta_hat(1); 1];                   % true mean tastes on constant, X1, dummy for segment
nrc           = size(theta_hat,1);
prods = size(bank_char_hat,1);
nlin          = 1;
theta   = theta_hat;
nmkt          = 10;
nobs          = nmkt*prods;                       % number of observations


% qv=[0:.1:1];
% qv=qv-mean(qv); 
% Nnodes=length(qv);
% qweight = repmat(1/Nnodes,Nnodes,1)';
% qweight = repmat(qweight,nobs,1);


% qv = 0:.1:1;
% qweight = pdf(makedist('Triangular','a',0,'b',1,'c',1),qv);
% qweight = qweight/sum(qweight);
% qv=qv-sum(qv.*qweight);
% Nnodes=length(qv);

nodes=length(qv);


id_g = bank_char_hat(:,3);
id_g          = repmat(id_g,nmkt,1);
id_t          = kron([1:nmkt]',ones(prods,1));

id_tg = grp2idx(num2str([id_t id_g]));
dum_tg = sparse(dummyvar(id_tg));
dum_t = sparse(dummyvar(id_t));


X1            = repmat(bank_char_hat(:,1),nmkt,1);
mc            = repmat(bank_char_hat(:,2),nmkt,1);



%  Data Structure
ffrgrid=linspace(0,20,nmkt)'*.5; % (min,max,nodes)
ffr=kron(ffrgrid,ones(prods,1));
Data.ffr            = ffr;
Data.nmkt           = nmkt;
Data.nobs           = nobs;
Data.theta    = theta;
Data.nodes         = nodes;
Data.nlin           = nlin;
Data.nrc            = nrc;
Data.X1             = X1;
Data.mc             = mc;
Data.betatrue       = betatrue;
Data.id_g           = id_g;
Data.qweight          = qweight;
Data.id_tg     = id_tg;
Data.dum_tg     = dum_tg;
Data.id_t     = id_t;
Data.dum_t     = dum_t;
Data.qv = qv ;

x=ffr(id_g==2);

%%
fun= @(x)foc(x,Data);
options = optimoptions('fsolve','Display','iter','MaxFunctionEvaluations',100000);
[price_hat,fval,exitflag] = fsolve(fun,x,options);

price            = zeros(size(ffr));
price(id_g==1)   = ffr(id_g==1); % cash
price(id_g==2)   = price_hat; % dep
price(price>ffr)=ffr(price>ffr); % ZLB
Xexo            = [price X1];
% Xrandom=[price trans_dum];
Xrandom=[price ];

xv=zeros(nobs,nodes*nrc);
for i = 1:nrc
    xv(:,(i-1)*nodes+1:i*nodes)=bsxfun(@times,Xrandom(:,i),qv);
end
delta           = Xexo*betatrue;
% Pack new rate vector to calculate market share
Data.xv             = xv;
[~,sj] = ShareCalculation(theta,delta,Data);


Data.sj = sj;
Data.price = price;
Data.logobsshare = log(sj);
data_simu = Data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pricegrid2 = grpstats(price(id_g==2),id_t(id_g==2),@(x)mean(x,1));
sjgrid2 = grpstats(sj(id_g==2),id_t(id_g==2),@(x)sum(x,1));
h=figure(1);
hold on;
yyaxis left
plot(ffrgrid,pricegrid2,'-')
ylim([0 1.5*max(pricegrid2)])
yyaxis right
plot(ffrgrid,sjgrid2,'d-')
legend('cb spread','cb dep');

title([ ' \beta:' num2str(betatrue(1),2) ' \sigma_1:' num2str(theta(1),2) ...
    ' \sigma_2:' num2str(theta(end),2)  ...
    ' prods:' num2str(prods)])
hold off;


pricegrid1 = grpstats(price(id_g==1),id_t(id_g==1),@(x)mean(x,1));
sjgrid1 = grpstats(sj(id_g==1),id_t(id_g==1),@(x)sum(x,1));
h=figure(2);
hold on;
yyaxis left
plot(ffrgrid,pricegrid1,'-')
ylim([0 1.5*max(pricegrid1)])
yyaxis right
plot(ffrgrid,sjgrid1,'d-')
legend('ffr','cash');

title([ ' \beta:' num2str(betatrue(1),2) ' \sigma_1:' num2str(theta(1),2) ...
    ' \sigma_2:' num2str(theta(end),2)  ...
    ' prods:' num2str(prods)])
hold off;
