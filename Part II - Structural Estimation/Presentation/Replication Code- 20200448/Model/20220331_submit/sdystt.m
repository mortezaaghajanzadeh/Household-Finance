function [Dist_ELfA1] = sdystt(Par,Policy,isValid, Dist_ELfA0, Transition_fAEL)

%         Dist_A = ones(1,Par.nA)/Par.nA*Par.A_Trans^100;
%         Dist_f = ones(1,Par.nf)/Par.nf*Par.f_Trans^100;
% 
%        Transition_A  = kron(eye(Par.nE*Par.nL),Par.A_Trans); 
%        Dist_ELfA1 = [];
%        
%        for p = 1:Par.nf
%            Transition_EL = Transition_A * 0;
%            for q = 1:Par.nA             
%                index0 = (p-1)*Par.nA + q;
%                state_beg = ((1:Par.nE*Par.nL)-1)*Par.nA+q;
%                state_end =  (Policy(index0,:)-1)*Par.nA+q;
%                Transition_EL((state_end-1)*Par.nE*Par.nL*Par.nA+state_beg) = 1;
%            end
%            
%            
%            Transition_valid = reshape(isValid((p-1)*Par.nA+1:p*Par.nA,:),[],1); 
%            Transition_valid = repmat(Transition_valid,1,Par.nE*Par.nL*Par.nA);
%            
%            Transition_new = repmat(kron(ones(Par.nE*Par.nL,1)/Par.nE/Par.nL,Dist_A'),1,Par.nE*Par.nL*Par.nA);
%            Transition_EL = Transition_EL.*(Transition_valid) + Transition_new.*(1-Transition_valid);
%            
%            Transition_ELA = sparse(Transition_EL*Transition_A);
%            Dist_ELA = 1/(Par.nE*Par.nL*Par.nA)*ones(1,Par.nE*Par.nL*Par.nA)*Transition_ELA^10;
%            Dist_ELfA1 = [Dist_ELfA1;Dist_f(p)*reshape(Dist_ELA,Par.nA,[])];
%        end

       
       Transition_EL = zeros(Par.nE*Par.nL*Par.nf*Par.nA, Par.nE*Par.nL*Par.nf*Par.nA);
       fA_distribution = (ones(1,Par.nA*Par.nf)/(Par.nA*Par.nf))*Par.fA_Trans^1000;
 
       for q = 1:Par.nA*Par.nf
            state_beg = ((1:Par.nE*Par.nL)-1)*Par.nA*Par.nf+q;
            state_end =  (Policy(q,:)-1)*Par.nA*Par.nf+q;
            Transition_EL((state_end-1)*Par.nE*Par.nL*Par.nA*Par.nf+state_beg) = 1;
       end
       Transition_EL = sparse(Transition_EL);

       Dist_ELfA_temp = sparse(Dist_ELfA0);
       Dist_ELfA1 = (reshape(isValid,1,[]).*Dist_ELfA_temp)*Transition_EL*Transition_fAEL;
       Dist_ELfA1 = Dist_ELfA1/sum(Dist_ELfA1);
           
       distribution_iter = 0;
       while sum(sum(abs(Dist_ELfA_temp-Dist_ELfA1))) > 0.0001
           Dist_ELfA_temp = Dist_ELfA1;

           Dist_ELfA1 = (reshape(isValid,1,[]).*Dist_ELfA_temp)*Transition_EL*Transition_fAEL;
           Dist_ELfA1 = reshape(Dist_ELfA1,Par.nA*Par.nf,[]);
           Dist_ELfA1 = Dist_ELfA1./repmat(sum(Dist_ELfA1,2),1,Par.nE*Par.nL);
           Dist_ELfA1 = Dist_ELfA1.*repmat(fA_distribution',1,Par.nE*Par.nL);
           Dist_ELfA1 = reshape(Dist_ELfA1,1,[]);
           Dist_ELfA1 = Dist_ELfA1/sum(Dist_ELfA1);
            
           distribution_iter= distribution_iter+1;
           str   = [' dist_iter = ', num2str(distribution_iter), ' , error = ', num2str(sum(abs(Dist_ELfA_temp-Dist_ELfA1)))];
%            disp(str);
           if distribution_iter > 2000
               disp('Distribution does not converge')
               break;
           end
       end

        Dist_ELfA1 = reshape(Dist_ELfA1,1,[]);
        
end