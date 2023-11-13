function [rd_bar, r_bar, D_avg, B_avg] = forecast_panel(Par, solution)

    load('shocks_6_500000_1_0_0'); 
    [t,N] = size(shock0.indexA);
    [firm, ~, ~] = SimulatePanel_short(Par, solution, shock0);
       
     rd_bar1= zeros(t, Par.nq);
     mean0  = sum(Par.qv.*Par.qweight);
     sensitivity0 = Par.alphaD + Par.sigalphaD*(Par.qv - mean0); sensitivity0 = sensitivity0 + (1-sensitivity0).*(sensitivity0<1);

     for j = 1:t
         exp_rate = [];
         log_rate = [];
         for i = 1:Par.nq        
                temp = exp(sensitivity0(i)*firm.rd1(j,:) + Par.lD);
                exp_rate(i) = ones(1,N)/N/sum(ones(1,N)/N)*temp'; 
                log_rate(i)= (log(exp_rate(i)) - Par.lD)/(sensitivity0(i));

                if sensitivity0(i) < 0  log_rate(i) = Par.rd_star; end
          end
             rd_bar1(j,:)   = log_rate;     
     end

     rd_bar0 = solution.rd_bar(:,shock0.indexf); rd_bar0 = rd_bar0';
%      [b,bint,r,rint,stats] = regress(reshape(rd_bar0,[],1),reshape(rd_bar1,[],1));
%      str = [' Coff = ', num2str(b), ' , R-sqr = ', num2str(stats(1))]; disp(str)

     stats = 1-(sqrt(sum(rd_bar1-rd_bar0).^2)/sqrt(sum(rd_bar1).^2));
     str = [' Coff = ', num2str(1), ' , R-sqr = ', num2str(stats(1))]; disp(str)

     r_bar1= zeros(t, 1);
     for j = 1:t        
          temp = exp(Par.nuB - Par.alpha*reshape(firm.r1(j,:),1,[]));
          exp_rate = ones(1,N)/N/sum(ones(1,N)/N)*temp';  
          log_rate = (log(exp_rate) - Par.nuB)/(- Par.alpha);
          r_bar1(j,:) = log_rate;     
     end

     r_bar0 = solution.r_bar(:,shock0.indexf); r_bar0 = r_bar0';
%      [b,bint,r,rint,stats] = regress(reshape(r_bar0,[],1),reshape(r_bar1,[],1));
%      str = [' Coff = ', num2str(b), ' , R-sqr = ', num2str(stats(1))]; disp(str)
            
     stats = 1-(sqrt(sum(r_bar1-r_bar0).^2)/sqrt(sum(r_bar1).^2));
     str = [' Coff = ', num2str(1), ' , R-sqr = ', num2str(stats(1))]; disp(str)

     rd_bar  = zeros(Par.nq,Par.nf);
     r_bar    = zeros(1,Par.nf);
     B_avg   = zeros(1,Par.nf);
     D_avg   = zeros(1,Par.nf);
       
     for i = 1:solution.nf
         indexi = find(firm.indexf0(:,1) == i);
         if isempty(indexi)
             rd_bar(:,i) = repmat(Par.rd_star(i),Par.nq,1);
             r_bar(i)  = Par.r_star(i);
             B_avg(i)  = Par.B_star(i);
             D_avg(i)  = Par.D_star(i);
         else
             rd_bar(:,i) = reshape(mean(rd_bar1(indexi,:)),[],1);
             r_bar(i)  = mean(r_bar1(indexi));
             B_avg(i)  = mean(mean(firm.B1(indexi,:)));
             D_avg(i)  = mean( mean(firm.D1(indexi,:)));
         end
     end
end