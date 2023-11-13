function [shock] = genshocks(Par, n, m, N, t)


% This function generates shocks for subsequent simulation
% The function can have two versions: if flag == 0, then all shocks are
% generated from the law of motion specified in the paper starting from t=1

% The function can also be used to generate impulse response, in which case
% the shocks equals n for the first N/2 periods, and it equals m for the
% second half

    if nargin <2
        m = 0;
        n = 0;
    end
    
    if nargin < 4
        load('Data\saved_seed');
        [t,N] = size(saved_seed.n1);
        randf  = saved_seed.u1(:,1);
        randA = saved_seed.u2;
        indexE = min(floor( saved_seed.u3*Par.nE)+1,Par.nE);
        indexL = min(floor( saved_seed.u4*Par.nL)+1,Par.nL);
    else
        randf = rand(t,1);
        randA = rand(t,N);
        indexE = randi([1 Par.nE],t,N);
        indexL = randi([1 Par.nL],t,N);
    end


    indexf = ones(t,1);
    indexA = ones(t,N);
    
    f_benchmark = cumsum(Par.f_distr);
    A_benchmark = cumsum(Par.A_distr);
    
    indexf(1) = indexf(1) + sum(randf(1) > f_benchmark);
    
    for Atemp = 1:Par.nA-1
        indexA(1,:) = indexA(1,:) + (randA(1,:) > A_benchmark(Atemp));
    end

    for i = 2:t
        f_benchmark = Par.f_Trans(indexf(i-1),:); f_benchmark = cumsum(f_benchmark);
        A_benchmark = Par.A_Trans(indexA(i-1,:),:); A_benchmark = A_benchmark'; A_benchmark = cumsum(A_benchmark);
        
        indexf(i) = indexf(i) + sum(randf(i) > f_benchmark);
    
        for Atemp = 1:Par.nA-1
            indexA(i,:) = indexA(i,:) + (randA(i,:) > A_benchmark(Atemp,:));
        end        
    end

    
   if m > 0 && n > 0
      tempf = (Par.fGrid(1:end-1) + Par.fGrid(2:end))/2;
      indexf(1:floor(t/2)) = sum(n>tempf)+1;
      indexf(floor(t/2)+1:end) = sum(m>tempf)+1;
   end
   
   tempf = Par.f_Trans*Par.fGrid'; tempf = tempf';
   deltaf = Par.fGrid(indexf(2:end)) - tempf(indexf(1:end-1));
   shock.deltaf = [0;deltaf'];
   
    shock.randA  = randA;
    shock.randf  = randf;
    shock.indexA = indexA;
    shock.indexf = indexf;
    
    
    shock.valueA = Par.AGrid(indexA);
    shock.valuef = Par.fGrid(indexf); shock.valuef=shock.valuef';
    shock.indexE = indexE;
    shock.indexL = indexL;    
    shock.valueE = Par.EGrid(indexE);
    shock.valueL = Par.LGrid(indexL);
    
    
    
end

