function [x,w] = ghqnorm(n,mu,s2)
% ghqnorm  Gauss-Hermite quadrature for the Normal distribution.
%
% Args:
%   n: Number of nodes
%   mu: Mean
%   s2: Variance

    [x0,w0] = gausshermite(n);
    x = x0*sqrt(2*s2) + mu;
    w = w0 / sqrt(pi);
end


function [x,w] = gausshermite(n)
% Gauss Hermite nodes and weights following "Numerical Recipes for C"
    MAXIT = 10;
    EPS   = 3e-14;
    PIM4  = 0.7511255444649425;
    
    x = zeros(n,1);
    w = zeros(n,1);
    
    m = floor(n+1)/2;
    for i=1:m
        if i == 1
            z = sqrt((2*n+1)-1.85575*(2*n+1)^(-0.16667));
        elseif i == 2
            z = z - 1.14*(n^0.426)/z;
        elseif i == 3
            z = 1.86*z - 0.86*x(1);
        elseif i == 4
            z = 1.91*z - 0.91*x(2);
        else
            z = 2*z - x(i-2);
        end
        for iter = 1:MAXIT
            p1 = PIM4;
            p2 = 0;
            for j=1:n
                p3 = p2;
                p2 = p1;
                p1 = z*sqrt(2/j)*p2 - sqrt((j-1)/j)*p3;
            end
            pp = sqrt(2*n)*p2;
            z1 = z;
            z = z1 - p1/pp;
            if abs(z-z1) <= EPS, break, end
        end
        if iter>MAXIT, error('too many iterations'), end
        x(i)     = z;
        x(n+1-i) = -z;
        w(i)     = 2/pp/pp;
        w(n+1-i) = w(i);
    end
    x(:) = x(end:-1:1);
end