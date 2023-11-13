clear;
clc;
filepath = '..\';
addpath(strcat(filepath,'blp\blp_subroutin'));
texfolder=strcat(filepath,'output\Tables');
outputfolder=strcat(filepath,'\blp\intermediate_output\');



for sample_num = 1:5
%%

samplelist={'dep_1','dep_2','dep_3','dep_4', 'dep_5'}; 
data = csvread(strcat('blp_matlab_data\',samplelist{sample_num},'.csv'),1);

id_t = data(:,1);
id_j = data(:,2); 
id_g = data(:,3);
id_tg = grp2idx(num2str([id_t id_g]));
dum_g = sparse(dummyvar(id_g)); 
dum_j = sparse(dummyvar(id_j)); 
dum_t = sparse(dummyvar(id_t));
dum_tg = sparse(dummyvar(id_tg));

A0 = [data(:,7:8)]; % attributes
price = data(:,5);
price(price<0)=0; % adding this filter helps to get rid of the spikes after rate decrease (does not change estimates much when using constrained estimation)
share = data(:,4);
logobsshare=log(share);
% z = data(:,9:12); % with shadow bank
z = data(:,9:10); % wo shaodw bank
nobs = size(id_t,1);
nmkt = size(unique(id_t),1) ;     % number of markets 
nbrn = size(unique(id_j),1) ; % number of brands
nseg = size(unique(id_g),1) ; % number of industry


% Matrices for Estimation
A=sparse([A0 dum_j(:,1:end) dum_t(:,2:end) ]); %need to remove year==1

Xexo=[price A];
Xrandom=[price ];
nlin=size(Xexo,2);
nrc=size(Xrandom,2);

[qv qweight] = nwspgr('KPU', 1, 10);
nodes=length(qv);
qweight = qweight';
qv = -qv';

xv=zeros(nobs,nodes*nrc);
for i = 1:nrc
    xv(:,(i-1)*nodes+1:i*nodes)=bsxfun(@times,Xrandom(:,i),qv);
end

Z=[A z z.^2 ];

% Weighting Matrix
norm=mean(mean(Z'*Z),2);
W=full(inv((Z'*Z)/norm)/norm);
xzwz            = Xexo'*Z*W*Z';
xzwzx           = xzwz * Xexo;
invxzwzx        = xzwzx\eye(size(xzwzx,2));



% Pack
Data.logobsshare             = logobsshare;
Data.xzwz       = xzwz;
Data.invxzwzx   = inv(xzwzx);
Data.Z          = Z;
Data.W          = W;
Data.Xexo       = Xexo;
Data.xv             = xv;
Data.nmkt           = nmkt;
Data.nobs           = nobs;
Data.nodes         = nodes;
Data.nlin           = nlin;
Data.nrc            = nrc;
Data.qweight          = qweight;
Data.id_tg     = id_tg;
Data.dum_tg     = dum_tg;
Data.id_t     = id_t;
Data.dum_t     = dum_t;
Data.qv = qv ;



%% Random Coefficient Nested Logit

t1      = cputime;

% OPTIMIZATION
delta0 = zeros(nobs,1);
save mvalold delta0

angmm = @(theta20)gmm_con(theta20,Data);

options = optimset( 'Display','iter',...
    'GradObj','on','TolCon',1E-6,...
    'TolFun',1E-6,'TolX',1E-6,...
    'Hessian', 'off','DerivativeCheck','off','FinDiffType','central');

theta20  = [ones(nrc,1)];     % starting values - just pick a low value of ;ambda to avoid overflow problems from the start
x_L     = [0*ones(nrc,1)];       % lower bound (only for the nested logit, to avoid numerical problems when pho gets negative)
x_U     = [10*ones(nrc,1)];       % upper bound (only for the nested logit, to avoid numerical problems when pho gets close to 1)
[theta_hat, FctVal, exitflag,~,~] = ...
    fmincon(angmm,theta20, [], [], [], [], x_L, x_U, [], options);


delta_hat = DeltaCalculation(theta_hat,Data);
beta_hat  = invxzwzx*(xzwz*delta_hat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8. Standard Errors of Demand Estimation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%BLP
delta_hat = DeltaCalculation(theta_hat,Data);
beta_hat  = invxzwzx*(xzwz*delta_hat);

xi_gmm = delta_hat-Xexo*beta_hat;
dgf=(size(Xexo,1)-size(Xexo,2));
ser_gmm=(xi_gmm'*xi_gmm)./dgf;
Rsq=1-var(xi_gmm)/var(delta_hat);
Rsq_gmm=1-(1-Rsq)*(nobs-1)/(dgf-1);

[se12gmm varcovar]=seblp(beta_hat,theta_hat,Data);

% calculating se of sum q
Dq_beta2=sum(exp(Xexo(id_g==2,2:end)*beta_hat(2:end)).*Xexo(id_g==2,2:end))./sum(exp(Xexo(id_g==2,2:end)*beta_hat(2:end)));
varq2=Dq_beta2*varcovar(2:end-1,2:end-1)*Dq_beta2';
seq2=sqrt(diag(varq2));

Dq_beta1=sum(exp(Xexo(id_g==1,2:end)*beta_hat(2:end)).*Xexo(id_g==1,2:end))./sum(exp(Xexo(id_g==1,2:end)*beta_hat(2:end)));
varq1=Dq_beta1*varcovar(2:end-1,2:end-1)*Dq_beta1';
seq1=sqrt(diag(varq1));
clear varcovar 


%OLS
ou_temp=1-dum_t'*share;
ou=ou_temp(id_t,:);
y=log(share)-log(ou);
bols=(Xexo'*Xexo)\(Xexo'*y);
est=y-Xexo*bols;
dgf=(size(Xexo,1)-size(Xexo,2));
Rsq=1-var(est)/var(y);
Rsq_ols=1-(1-Rsq)*(nobs-1)/(dgf-1);

varcovar=est'*est/dgf*inv(Xexo'*Xexo);
seols=sqrt(diag(varcovar));
clear varcovar 


%TSLS
btsls=invxzwzx*(Xexo'*Z*W*Z'*y);
xi=y-Xexo*btsls;
Rsq=1-var(xi)/var(y);
Rsq_tsls=1-(1-Rsq)*(nobs-1)/(dgf-1);
varcovar=xi'*xi/nobs*invxzwzx;
setsls=sqrt(diag(varcovar));

clear varcovar 


%%

% implied semi-elasticity
[OwnSemiElast_hat] = OwnElast(beta_hat(1),theta_hat,delta_hat,Data);

sum_share_tg = accumarray(id_tg,share, [], @sum);
sum_share_tg=sum_share_tg(id_tg,:);
weight = share./sum_share_tg;
id_g_tg = accumarray(id_tg,id_g, [], @mean); % index
 

num_rep_bank = 6;
bank_char_hat=zeros(num_rep_bank+1,3); % quality, mc, industry
sum_exp_delta = accumarray(id_tg,exp(delta_hat), [], @sum); % sum 
sum_delta = accumarray(id_g_tg,log(sum_exp_delta), [], @mean); % average across markets
mean_price_tg = accumarray(id_tg,price.*weight, [], @sum);
mean_price_g = accumarray(id_g_tg,mean_price_tg, [], @mean);
q_sum = sum_delta - mean_price_g.*(beta_hat(1)+mean(qv).*theta_hat(1));
q_bank = q_sum - log([1;num_rep_bank]);
bank_char_hat(:,1) =q_bank([1;repmat(2,num_rep_bank,1)],:);

mc_hat = price-(-1./OwnSemiElast_hat); 
mean_mc_hat_tg = accumarray(id_tg,mc_hat.*weight, [], @sum);
mean_mc_hat_g = accumarray(id_g_tg,mean_mc_hat_tg, [], @mean);
bank_char_hat(:,2) = mean_mc_hat_g([1;repmat(2,num_rep_bank,1)],:); % heter bank char 
bank_char_hat(:,3) =[1;repmat(2,num_rep_bank,1)];

exp(sum_delta(1))/exp(sum_delta(2)) % ratio of market share when price>0
exp(q_sum(1))/exp(q_sum(2)) % ratio of market share when price=0 


parameters={};
parameters{1,1}='alpha';parameters{1,2}=num2str(beta_hat(1)+sum(qv.*qweight)*theta_hat(1));
parameters{2,1}='sigma_alpha';parameters{2,2}=num2str(theta_hat(1)/sqrt(12));%divide by sqrt(12) to convert the range of a uniform distribution to standard deviation
parameters{3,1}='q_cash';parameters{3,2}=num2str(q_sum(1,1));
parameters{4,1}='q_deposit_sum';parameters{4,2}=num2str(q_sum(2,1));
parameters{5,1}='N';parameters{5,2}=num2str(num_rep_bank(1));
parameters{6,1}='semi-elasticity';parameters{6,2}=num2str(mean(OwnSemiElast_hat(id_g==2)));
parameters{7,1}='elasticity';parameters{7,2}=num2str(mean(OwnSemiElast_hat(id_g==2).*price(id_g==2)));
parameters{8,1}='marginal cost';parameters{8,2}=num2str(mean_mc_hat_g(2));
parameters{9,1}='alpha (se)';parameters{9,2}=num2str(se12gmm(1));
parameters{10,1}='sigma_alpha (se)';parameters{10,2}=num2str(se12gmm(end)/sqrt(12));%divide by sqrt(12) to convert the range of a uniform distribution to standard deviation
parameters{11,1}='q_cash (se)';parameters{11,2}=num2str(seq1);
parameters{12,1}='q_deposit_sum (se)';parameters{12,2}=num2str(seq2);

output_table_name=strcat(outputfolder,'\demand_',samplelist{sample_num},'.xls');
xlswrite(output_table_name, parameters)

nchar=size([price A0],2);
save(samplelist{sample_num},'beta_hat','theta_hat','nchar','nrc','se12gmm', 'Rsq_gmm','nobs','qv','qweight') % the deposit file compile both results in one table


end

