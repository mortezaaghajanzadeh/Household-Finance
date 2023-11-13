function [ID, infl, count] = getinfl(filename)


% mmt0= [simu_un0.meanC2V, ...
%              simu_un0.meanN2D, ...
%              simu_un0.stdN2D, ...
%              simu_un0.spreadrd, ...
%              simu_un0.spreadr, ...
%              simu_un0.meanD2A, ...
%              simu_un0.meanFIX1 - Par.meanA, ...
%              simu_un0.meanLEV, ...
%              simu_un0.meanM2B, ...
%              simu_un0.reg_T2T_agg*(1-Par.mu)];
         
data = readtable(filename);
Nob  = size(data,1);

ID = data(:,2);  ID = table2array(ID); 

weight = data(:,12);  weight = table2array(weight); 
sum_weight = sum(weight);


mC2V = zeros(Nob,1);
datai = data(:,3);  datai = table2array(datai).*weight*Nob/sum_weight ;
indexi = find((datai < 1e4).*(datai > -1e4));
datai_indexi = datai(indexi) - mean(datai(indexi)); 
mC2V(indexi) = datai_indexi; count(1) = length(datai_indexi);

mN2D = zeros(Nob,1);
datai = data(:,5);  datai = table2array(datai).*weight*Nob/sum_weight ;
indexi = find((datai < 1e4).*(datai > -1e4));  mean_weighted_N2D = mean(datai(indexi));
datai_indexi = datai(indexi) - mean_weighted_N2D;  
mN2D(indexi) = datai_indexi; count(2) = length(datai_indexi);

sN2D = zeros(Nob,1);
datai = data(:,5);  datai = table2array(datai).^2.*weight*Nob/sum_weight - mean_weighted_N2D^2;
indexi = find((datai < 1e4).*(datai > -1e4));  
datai_indexi = datai(indexi) - mean(datai(indexi)); 
sN2D(indexi) = datai_indexi; count(3) = length(datai_indexi);
sN2D(indexi) = sN2D(indexi)/2/sqrt(mean(datai(indexi)));

mspreadrd = zeros(Nob,1);
datai = data(:,6);  datai = table2array(datai)/100.*weight*Nob/sum_weight;
indexi = find((datai < 1e4).*(datai > -1e4));  
datai_indexi = datai(indexi) - mean(datai(indexi)); 
mspreadrd(indexi) = datai_indexi; count(4) = length(datai_indexi);

mspreadr = zeros(Nob,1);
datai = data(:,7);  datai = table2array(datai)/100.*weight*Nob/sum_weight;
indexi = find((datai < 1e4).*(datai > -1e4));  
datai_indexi = datai(indexi) - mean(datai(indexi)); 
mspreadr(indexi) = datai_indexi; count(5) = length(datai_indexi);

mD2A = zeros(Nob,1);
datai = data(:,8);   datai = table2array(datai).*weight*Nob/sum_weight ;
indexi = find((datai < 1e4).*(datai > -1e4));  
datai_indexi = datai(indexi) - mean(datai(indexi)); 
mD2A(indexi) = datai_indexi; count(6) = length(datai_indexi);

mFIX = zeros(Nob,1);
datai = data(:,9);   datai = table2array(datai).*weight*Nob/sum_weight ;
indexi = find((datai < 1e4).*(datai > -1e4));  
datai_indexi = datai(indexi) - mean(datai(indexi)); 
mFIX(indexi) = datai_indexi; count(7) = length(datai_indexi);

mLEV = zeros(Nob,1);
datai = data(:,10);   datai = table2array(datai).*weight*Nob/sum_weight ;
indexi = find((datai < 1e4).*(datai > -1e4));  
datai_indexi = datai(indexi) - mean(datai(indexi)); 
mLEV(indexi) = datai_indexi; count(8) = length(datai_indexi);

mM2B = zeros(Nob,1);
datai = data(:,11);    datai = table2array(datai).*weight*Nob/sum_weight ;
indexi = find((datai < 1e4).*(datai > -1e4));  
datai_indexi = datai(indexi) - mean(datai(indexi)); 
mM2B(indexi) = datai_indexi; count(9) = length(datai_indexi);

infl = [mC2V, mN2D, sN2D, mspreadrd, mspreadr, mD2A, mFIX, mLEV, mM2B];

end