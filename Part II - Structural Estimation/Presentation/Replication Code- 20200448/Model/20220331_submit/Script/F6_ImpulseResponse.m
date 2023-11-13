%% First, expand solutions becuase we need finer grids in the impulse response
clear all
% main
load('Results\main')  % can directly load solution; otherwise, run main.m

Par.mf         = 5;	%	
[Par_expand,solution_expand] = expand_solution(Par,solution0);

%% Now, simulate transitions from old to new steady state

N    = 1e5;
t    = 100;
flag = 1;

f_pre  = 0.0086;   % we choose 0.86% becuase it is the "flection point" in Figure 5, Panel A
f_post = 0.02;      % we examine the impulse reponse of a decrese in FFR. 

shock1= genshocks(Par_expand, f_pre, f_post, N, t); 
[banks1, simu_con1, simu_un1, rec1] = SimulatePanel_short(Par_expand, solution_expand, shock1);


N    = 1e5;
t    = 100;
flag = 1;

f_pre  = 0.0086;  % we choose 0.86% becuase it is the "flection point" in Figure 5, Panel A
f_post = 0.001;   % we examine the impulse reponse of a decrese in FFR. 

shock2= genshocks(Par_expand, f_pre, f_post, N, t); 
[banks2, simu_con2, simu_un2, rec2] = SimulatePanel_short(Par_expand, solution_expand, shock2);

save("Results\F6_ImpulseResponse",'banks1','banks2', '-v7.3')

%% Plot the Results 

load("Results\F6_ImpulseResponse")
for i  = [1,2]
eval(['banks  = banks',num2str(i),';']);

Dt = mean(banks.D1,2);
Bt  = mean(banks.B1,2);
Lt  = mean(banks.L1,2);
Nt  = mean(banks.N1,2);
E0t = mean(banks.E0,2);

Bank_Capt = E0t;

% use the last period pre-transition values to measure the "old" steady state
Dt(1:49) = Dt(50);
Bt(1:49) = Bt(50);
Lt(1:49) = Lt(50);
Bank_Capt(1:49) = Bank_Capt(50);


x =49:60;

figure(6)
subplot(2,1,i)
set(0, 'DefaultAxesFontSize', 12) 
set(0,'defaultTextInterpreter','latex');

[ax,p1,p2] = plotyy([x',x'],[Lt(x)/Lt(50),Dt(x)/Dt(50)],[x',x'],[Bank_Capt(x)/Bank_Capt(50),Nt(x)/Nt(50)]); 
xlim(ax(1),[49, 60])
xlim(ax(2),[49, 60])

set(p1(1),'color',[0.8500, 0.3250, 0.0980],'marker','o','MarkerSize',5,'MarkerFaceColor',[0.8500, 0.3250, 0.0980])
set(p1(2),'color',[0, 0.4470, 0.7410],'marker','d','MarkerSize',5,'MarkerFaceColor',[0, 0.4470, 0.7410])
set(p2(1),'color',[0.9290, 0.6940, 0.1250],'marker','*')
set(p2(2),'color',[0.4940, 0.1840, 0.5560],'marker','+')

if i ==1
    ylim(ax(1),[0.8 1.2])
    ylim(ax(2),[0.8 1.2])
    set(ax(1),'ytick',linspace(0.8, 1.2, 5));
    set(ax(2),'ytick',linspace(0.8, 1.2, 5));
    title("Panel A: Impulse response to increases in the Federal funds rate"+newline+"   ",'FontSize', 12,'FontWeight', 'Normal')
elseif i ==2
    ylim(ax(1),[0.85 1.05])
    ylim(ax(2),[0.4 1.2])
    set(ax(1),'ytick',linspace(0.85, 1.05, 5));
    set(ax(2),'ytick',linspace(0.4, 1.2, 5));
    title("Panel B: Impulse response to decreases in the Federal funds rate"+newline+"   ",'FontSize', 12,'FontWeight', 'Normal')

end

set(ax(1),'fontsize',12)
set(ax(2),'fontsize',12)
xticks([50 51 52 53 54 55 56 57 58 59 60])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
 
ylabel(ax(1), 'Bank Loans/Deposit ','FontSize', 12);
ylabel(ax(2), 'Capital/Non-reservables', 'FontSize', 12);

end

