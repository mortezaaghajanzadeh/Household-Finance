function [X,V,c] = solveLifeCycleHH(r,Nw,Nr,genPars,incPars,method)
    % SOLVELIFECYCLEHH  Solve the household's life-cycle model.
    %
    % Args:
    %   r: Interest rate on savings
    %   Nw: Number of working-life periods
    %   Nr: Number of retirement periods
    %   genPars = {a,beta,gamma,ps}: General model parameters, with
    %       a: Asset holdings grid
    %       beta: Household's discount factor
    %       gamma: Household's relative risk aversion coefficient
    %       ps: (Nw+Nr)-by-1 vector of probabilities of survival
    %   incPars = {G,yfloor,psi,yp,Pp,yt,Pt}: Income process parameters, with
    %       G: Nw-by-1 vector of life-cycle 'base' income levels
    %       yfloor: Floor on disposable income
    %       psi: Replacement rate for retirement income
    %       yp, Pp: Permanent income grid and corresponding transition matrix
    %       yt, Pt: Transitory income grid and corresponding probabilities
    %   method: Interpolation method for EGM
    %
    [a,beta,gamma,ps] = genPars{:};
    [G,yfloor,psi,yp,Pp,yt,Pt] = incPars{:};
    N = Nw + Nr;       % Total number of periods
    Na = size(a,1);    % Asset grid size
    amin = a(1);       % Borrowing constraint
    Nyp = size(yp,1);  % Number of permanent income states
    Nyt = size(yt,1);  % Number of transitory income states
    V = cell(1,N);     % Value function
    c = cell(1,N);     % Consumption policy
    
    % Create cash-on-hand grids
    X = cell(1,Nw);
    ls_yt = linspace(yt(1),yt(Nyt),Na)';
    for j = 1:Nw
        X{j} = (1+r)*a + max(G(j) * yp' .* ls_yt,yfloor);
    end
    for j = Nw+1:N
        % Keep track of permanent income state at j=Nw
        X{j} = (1+r)*a + psi * max(G(Nw) * yp',yfloor);
    end
    
    % A. Retirement period (Nw -> N) [No income risk]
    c{N} = X{N};
    V{N} = c{N}.^(1-gamma) / (1-gamma);
    for j = N-1:-1:Nw
        Va = (1+r) * c{j+1}.^(-gamma);
        cj = (beta*ps(j+1) * Va).^(-1/gamma);
        xj = cj + a;  % Endogenous cash-on-hand grid
        Vj = cj.^(1-gamma) / (1-gamma) + beta*ps(j+1) * V{j+1};
        for s = 1:Nyp
            xjk = xj(:,s);
            cjk = cj(:,s);
            if amin < xj(1,s)  % Add x = a(1), c = 0 point to interpolation set
                xjk = [amin;xjk];
                cjk = [0;cjk];
            end
            cj(:,s) = interp1(xjk,cjk,X{j}(:,s),method,'extrap');
            Vj(:,s) = interp1(xj(:,s),Vj(:,s),X{j}(:,s),method,'extrap');
        end
        c{j} = max(min(cj,X{j}-amin),0);  % Enforce constraints
        V{j} = Vj;
    end

    % B. Working-life period (1 -> Nw-1)
    for j = Nw-1:-1:1
        EV1 = zeros(Na,Nyp);  % Expectation over transitory income shock
        EVa = zeros(Na,Nyp);
        for s = 1:Nyp
            V1p = zeros(Na,Nyt);
            c1p = zeros(Na,Nyt);  % c'(a,re) for y fixed
            for k = 1:Nyt
                Xpt = (1+r)*a + max(G(j+1)*yp(s)*yt(k),yfloor);
                V1p(:,k) = interp1(X{j+1}(:,s),V{j+1}(:,s),Xpt,method,'extrap');
                c1p(:,k) = interp1(X{j+1}(:,s),c{j+1}(:,s),Xpt,method,'extrap');
            end
            EV1(:,s) = V1p * Pt;
            EVa(:,s) = ((1+r) * c1p.^(-gamma)) * Pt;
        end
        cj = (beta*ps(j+1) * EVa * Pp').^(-1/gamma);
        xj = cj + a;  % Endogenous cash-on-hand grid
        Vj = cj.^(1-gamma) / (1-gamma) + beta*ps(j+1) * EV1 * Pp';
        for s = 1:Nyp
            xjk = xj(:,s);
            cjk = cj(:,s);
            if amin < xj(1,s)  % Add x = a(1), c = 0 point to interpolation set
                xjk = [amin;xjk];
                cjk = [0;cjk];
            end
            cj(:,s) = interp1(xjk,cjk,X{j}(:,s),method,'extrap');
            Vj(:,s) = interp1(xj(:,s),Vj(:,s),X{j}(:,s),method,'extrap');
        end
        c{j} = max(min(cj,X{j}-amin),0);  % Enforce borrowing constraint
        V{j} = Vj;
    end
end
