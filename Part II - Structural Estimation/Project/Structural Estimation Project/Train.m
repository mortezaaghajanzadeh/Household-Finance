%%

clear
close all
clc

useGPU = 0;

%% training set

nSimul = 100000;

Allp = haltonset(3,'Skip',1);  Allp = scramble(Allp,'RR2');  Allp = net(Allp,1000) ;
LB = [.8 1.01 0]; UB = [.99 5 .01];
Allp = LB + Allp.*(UB-LB);
save("data/Allp.mat","Allp")

%% train
for thiS =  1:size(Allp,1)

disp(thiS)
disp(Allp(thiS,:))

m = simul(Allp(thiS,:),nSimul,useGPU,0);



disp(m)

Allm(thiS,:) = m;

disp('------------------------')

end

save("data/Allm.mat","Allm")





