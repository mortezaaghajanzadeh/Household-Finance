function [Z,Zprob] = tauchenHussey(N,mu,rho,sigma,baseSigma)
% Function TAUCHENHUSSEY
%
% Purpose:    Finds a Markov chain whose sample paths
%             approximate those of the AR(1) process
%                 z(t+1) = (1-rho)*mu + rho * z(t) + eps(t+1)
%             where eps are normal with stddev sigma
%
% Format:     {Z, Zprob} = TauchenHussey(N,mu,rho,sigma,m)
%
% Input:      N         scalar, number of nodes for Z
%             mu        scalar, unconditional mean of process
%             rho       scalar
%             sigma     scalar, std. dev. of epsilons
%             baseSigma scalar, std. dev. used to calculate Gaussian
%                       quadrature weights and nodes, i.e. to build the
%                       grid. I recommend that you use baseSigma = w*sigma +
%                       (1-w)*sigmaZ where sigmaZ = sigma/sqrt(1-rho^2),
%                       and w = 0.5 + rho/4. Tauchen & Hussey recommend
%                       baseSigma = sigma, and also mention baseSigma = sigmaZ.
%
% Output:     Z       N*1 vector, nodes for Z
%             Zprob   N*N matrix, transition probabilities
%
%     Martin Floden, Stockholm School of Economics
%     January 2007 (updated August 2007)
%
%     This procedure is an implementation of Tauchen and Hussey's
%     algorithm, Econometrica (1991, Vol. 59(2), pp. 371-396)

    Zprob = zeros(N,N);
    [Z,w] = ghqnorm(N,mu,baseSigma^2);   % See note 1 below
    
    for i = 1:N
        for j = 1:N
            EZprime    = (1-rho)*mu + rho*Z(i);
            Zprob(i,j) = (w(j) * normpdf(Z(j),EZprime,sigma) ...
                          / normpdf(Z(j),mu,baseSigma) ...
                          );
        end
    end
    for i = 1:N
        Zprob(i,:) = Zprob(i,:) / sum(Zprob(i,:),2);
    end
end

% Note 1: If you have Miranda and Fackler's CompEcon toolbox you can use
% their qnwnorm function to obtain quadrature nodes and weights for the
% normal function: [Z,w] = qnwnorm(N,mu,baseSigma^2);
% Compecon is available at http://www4.ncsu.edu/~pfackler/compecon/
% Otherwise, use gaussnorm as here. 
