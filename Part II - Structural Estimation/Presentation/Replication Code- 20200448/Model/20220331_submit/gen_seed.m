function [] = gen_seed(N,t)

% Generate and store random seed for simulation
% This function only needs to run once
% Seeds are kept constant throughout the life of the project

saved_seed.N = N;
saved_seed.t  = t;

t = t+2;

saved_seed.u1 = rand(t,N);
saved_seed.u2 = rand(t,N);
saved_seed.u3 = rand(t,N);
saved_seed.u4 = rand(t,N);
saved_seed.u5 = rand(t,N);
saved_seed.u6 = rand(t,N);

saved_seed.n1 = randn(t,N);
saved_seed.n2 = randn(t,N);
saved_seed.n3 = randn(t,N);
saved_seed.n4 = randn(t,N);
saved_seed.n5 = randn(t,N);
saved_seed.n6 = randn(t,N);

save('Data\saved_seed','saved_seed');

end