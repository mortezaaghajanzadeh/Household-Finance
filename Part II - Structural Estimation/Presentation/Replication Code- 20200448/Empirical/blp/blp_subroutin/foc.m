function f = foc(x,Data)
% Unpack
nmkt            = Data.nmkt;
theta     = Data.theta;
nodes          = Data.nodes;
qv          = Data.qv;
qweight   = Data.qweight;
betatrue        = Data.betatrue;
ffr             = Data.ffr;
nrc             = Data.nrc;
id_g            = Data.id_g;
X1              = Data.X1;
mc              = Data.mc;
nobs = Data.nobs;
% price vector
price            = zeros(size(ffr));
price(id_g==1)   = ffr(id_g==1); % cash
price(id_g~=1)   = x; % dep
trans_dum = (id_g<=2);

Xexo            = [price X1];
% Xrandom = [price trans_dum];
Xrandom = [price ];

xv=zeros(nobs,nodes*nrc);
for i = 1:nrc
    xv(:,(i-1)*nodes+1:i*nodes)=bsxfun(@times,Xrandom(:,i),qv);
end


alpha           = betatrue(1);
delta           = Xexo*betatrue;

% Pack new rate vector to calculate markup
Data.xv             = xv;


% Random Coefficient Nested Logit wrt rate
[OwnSemiElast] = OwnElast(alpha,theta,delta,Data);

%% foc
% f = ffr(id_s==2)-x-1./OwnSemiElast(id_s==2); % need to remove cash from FOC since FFR is exogenous
f = x+1./OwnSemiElast(id_g==2)-mc(id_g==2); % need to remove cash from FOC since FFR is exogenous