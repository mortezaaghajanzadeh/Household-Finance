function f=jacob(theta,delta,Data)

% jacob
% Jacobian of the implicit function that defines the mean utility
% This function is needed to compute the optimal instruments at a given
% value of theta
% 
% Written by Mathias Reynaert (2013)
% Original Source: Aviv Nevo (2000)

id_t = Data.id_t;
nrc = Data.nrc;
nodes = Data.nodes;
qweight = Data.qweight;
xv = Data.xv;

%% Market Shares
[sij,~]=ShareCalculation(theta,delta,Data);
wsij = qweight.*sij;

%% Jacobian
derTheta=zeros(size(xv,1),size(theta,1));
for t = 1:max(id_t)
wssij=wsij((id_t==t),:);
ssij=sij((id_t==t),:);
% sqweight=qweight((id_t==j),:);
part1=wssij*ssij';
derShareDelt = (diag(sum(wssij,2)) - part1);
f1 = zeros(size(derShareDelt,1),nrc);
% computing (partial share)/(partial sigma)
for i = 1:nrc   
 	sxv=xv(id_t==t,((i-1)*nodes)+1:i*nodes);
    sumxv=sum(sxv.*ssij,1);
    sumxv=repmat(sumxv,size(sxv,1),1);
%  	f1(:,i) = sum(sqweight.*(ssij.*(sxv-sumxv)),2);
 	f1(:,i) = sum(qweight.*(ssij.*(sxv-sumxv)),2);
 	clear sxv sumxv
end
derTheta((id_t==t),:)=-derShareDelt\f1;
end
f=derTheta;

end