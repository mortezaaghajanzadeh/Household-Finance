clear all
% main
load('Results\main')  % can directly load solution; otherwise, run main.m
load ('spread_data.mat') 

log_FFR = log(data(:,4));

 mean_ffr = mean(log_FFR);  
 std_ffr = std(log_FFR) ;         
 beta_ffr= regress(log_FFR(2:end)-mean_ffr,[log_FFR(1:end-1)-mean_ffr,log_FFR(1:end-1)*0+1]); 

Par_temp = Par;
Par_temp.rhof       = beta_ffr(1);                                                 	
Par_temp.nuf        = mean_ffr - log(100);                                          	
Par_temp.sigf       = sqrt(std_ffr*(1-Par.rhof^2));                                 	

[~,~, Par_temp.f_Trans] = fATransition(Par_temp);
  
  Ef = Par_temp.fGrid';
  for i = 1:100
      Ef  =  Ef + (1 - Par_temp.mu)^i*(Par_temp.f_Trans^i)*Par_temp.fGrid';
  end
 
  Eyrs = 1;
  for i = 1:100
      Eyrs  =  Eyrs + (1 - Par_temp.mu)^i*1;
  end

  EF5 = Ef/Eyrs;  % calculate the 5-year average FFR


loan_spread = data(:,2)/100 +Par.meanA;
deposit_spread = max(data(:,3)/100,0);
FFR = data(:,4)/100;
FFR5 = interp1(Par.fGrid,EF5,FFR,'linear','extrap');

loan_rate = loan_spread + FFR5;
deposit_rate = FFR - deposit_spread;
save("Results\F4_Compare_spreads",'FFR','deposit_rate','loan_rate','Par','simu_con0')

%% Plot results

load("Results\F4_Compare_spreads")
figure(4)
subplot(1,2,1)
r3 = ksrlin(FFR,deposit_rate);
scatter(FFR,deposit_rate,'MarkerEdgeColor',[0.5 0.5 0.5]); hold on
plot(r3.x, r3.f), hold on
plot(Par.fGrid, smooth(simu_con0.rd1),'b--'); hold on
xlim([0,0.08])
ylim([0,0.08])
xlabel("Federal funds rate"+newline+ "   "+newline+ "   ")
ylabel('Deposit rate')

subplot(1,2,2)
r4 = ksrlin(FFR,loan_rate);
scatter(FFR,loan_rate,'MarkerEdgeColor',[0.5 0.5 0.5]); hold on
plot(r4.x,r4.f), hold on
plot(Par.fGrid, smooth(simu_con0.r1,'loess'),'b--')
xlim([0,0.08])
ylim([0,0.08])
xlabel("Federal funds rate"+newline+ "   "+newline+ "   ")
ylabel('Loan rate')
legend('Raw data','Local polynomial smooth plots using raw data','Model predictions','location','southoutside','Orientation','horizontal'); legend boxoff

