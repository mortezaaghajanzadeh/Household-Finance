function [out, out_A, out_f] = fATransition(Par)

% This function calculate the bivariant transition matrix (f,A)--->(f',A')

% |f'-Par.nuf|   |Par.rhof,  0      |   |f-Par.nuf|  |N(0,sigf)|
% |          | = |                  | * |         |+ |         |
% |A'-Par.nuA|   |Par.rhoAf, Par.rhA|   |A-Par.nuA|  |N(0,sigA)|

out = -1*ones(Par.nf*Par.nA,Par.nf*Par.nA);

log_fGrid = log(Par.fGrid);
log_AGrid = log(Par.AGrid);

if Par.nA ~=1 && Par.nf ~=1
    
for i0 = 1:Par.nf
    for j0 = 1:Par.nA
        for i1 = 1:Par.nf
            for j1 = 1:Par.nA
                
               index0 = (i0-1)*(Par.nA)+j0;
               index1 = (i1-1)*(Par.nA)+j1;
               
               deltaf = log_fGrid(i1) - Par.nuf - Par.rhof*(log_fGrid(i0) - Par.nuf);
               deltaA = log_AGrid(j1) - Par.nuA - Par.rhoA*(log_AGrid(j0) - Par.nuA);
                                           
               if j1 == 1 
                   upperA = (log_AGrid(j1) + log_AGrid(j1+1))/2 -  Par.rhoA*(log_AGrid(j0) - Par.nuA) - Par.nuA;
                   ProbA  = normcdf(upperA,0,Par.sigA);
               elseif j1 == Par.nA
                   lowerA = (log_AGrid(j1) + log_AGrid(j1-1))/2 -  Par.rhoA*(log_AGrid(j0) - Par.nuA) - Par.nuA;
                   ProbA  = 1 - normcdf(lowerA,0,Par.sigA);
               else
                   upperA = (log_AGrid(j1) + log_AGrid(j1+1))/2 -  Par.rhoA*(log_AGrid(j0) - Par.nuA) - Par.nuA;
                   lowerA = (log_AGrid(j1) + log_AGrid(j1-1))/2 -  Par.rhoA*(log_AGrid(j0) - Par.nuA) - Par.nuA;
                   ProbA  = normcdf(upperA,0,Par.sigA) - normcdf(lowerA,0,Par.sigA);                  
               end
               
               if i1 == 1
                   upperf = (log_fGrid(i1) + log_fGrid(i1+1))/2 -  Par.rhof*(log_fGrid(i0) - Par.nuf) - Par.nuf;
                   Probf  = normcdf(upperf,0,Par.sigf);
               elseif i1 == Par.nf
                   lowerf = (log_fGrid(i1) + log_fGrid(i1-1))/2 -  Par.rhof*(log_fGrid(i0) - Par.nuf) - Par.nuf;
                   Probf  = 1 - normcdf(lowerf,0,Par.sigf);
               else
                   upperf = (log_fGrid(i1) + log_fGrid(i1+1))/2 -  Par.rhof*(log_fGrid(i0) - Par.nuf) - Par.nuf;
                   lowerf = (log_fGrid(i1) + log_fGrid(i1-1))/2 -  Par.rhof*(log_fGrid(i0) - Par.nuf) - Par.nuf;
                   Probf  = normcdf(upperf,0,Par.sigf) - normcdf(lowerf,0,Par.sigf);
               end
               
               out(index0,index1) = ProbA * Probf;            
               out_A(j0,j1) = ProbA;
               out_f(i0,i1) = Probf;
            end
        end
    end
end

elseif Par.nA ==1
    
   halff = (log_fGrid(2)-log_fGrid(1))/2;
   
    for i0 = 1:Par.nf
        for i1 = 1:Par.nf
                
               index0 = i0;
               index1 = i1;
               
               deltaf = log_fGrid(i1) - Par.nuf - Par.rhof*(log_fGrid(i0) - Par.nuf);               
               
               if i1 == 1 
                   upperf = (log_fGrid(i1) + log_fGrid(i1+1))/2 -  Par.rhof*(log_fGrid(i0) - Par.nuf) - Par.nuf;
                   Probf  = normcdf(upperf,0,Par.sigf);
               elseif i1 == Par.nf
                   lowerf = (log_fGrid(i1) + log_fGrid(i1-1))/2 -  Par.rhof*(log_fGrid(i0) - Par.nuf) - Par.nuf;
                   Probf  = 1 - normcdf(lowerf,0,Par.sigf);
               else
                   upperf = (log_fGrid(i1) + log_fGrid(i1+1))/2 -  Par.rhof*(log_fGrid(i0) - Par.nuf) - Par.nuf;
                   lowerf = (log_fGrid(i1) + log_fGrid(i1-1))/2 -  Par.rhof*(log_fGrid(i0) - Par.nuf) - Par.nuf;
                   Probf  = normcdf(upperf,0,Par.sigf) - normcdf(lowerf,0,Par.sigf);
               end
               
               out(index0,index1) = Probf;     
               out_f = out; out_A=1;
        end
    end
    
elseif Par.nf ==1
   
    halfA = (log_AGrid(2)-log_AGrid(1))/2;
    
    for j0 = 1:Par.nA
        for j1 = 1:Par.nA
                
               index0 = j0;
               index1 = j1;
               
               deltaA = log_AGrid(j1) - Par.nuA - Par.rhoA*(log_AGrid(j0) - Par.nuA);

               if j1 == 1 
                   upperA = (log_AGrid(j1) + log_AGrid(j1+1))/2 -  Par.rhoA*(log_AGrid(j0) - Par.nuA) - Par.nuA;
                   ProbA  = normcdf(upperA,0,Par.sigA);
               elseif j1 == Par.nA
                   lowerA = (log_AGrid(j1) + log_AGrid(j1-1))/2 -  Par.rhoA*(log_AGrid(j0) - Par.nuA) - Par.nuA;      
                   ProbA  = 1 - normcdf(lowerA,0,Par.sigA);
               else
                   upperA = (log_AGrid(j1) + log_AGrid(j1+1))/2 -  Par.rhoA*(log_AGrid(j0) - Par.nuA) - Par.nuA;
                   lowerA = (log_AGrid(j1) + log_AGrid(j1-1))/2 -  Par.rhoA*(log_AGrid(j0) - Par.nuA) - Par.nuA;
                   ProbA = normcdf(upperA,0,Par.sigA) - normcdf(lowerA,0,Par.sigA);                  
               end              
               
               out(index0,index1) = ProbA;
               out_A = out; out_f=1;       
        end
    end
    
end


end

