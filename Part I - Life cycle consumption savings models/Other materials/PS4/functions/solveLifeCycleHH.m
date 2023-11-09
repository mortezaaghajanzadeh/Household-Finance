function [V,c] = solveLifeCycleHH(Xw,Xr,P,pars,method)
    % SOLVELIFECYCLEHH  Solve the household's life-cycle model.
    %
    % Args:
    %   Xw, Xr: Cash-on-hand grids for working and retired households.
    %   P: Income transition matrix
    %   pars: Fixed model parameters = {a,beta,gamma,r,N,ps}
    %   method: Interpolation method for EGM
    %
    [a,beta,gamma,r,N,ps] = pars{:};
    Nw = size(Xw,2);
    amin = a(1);         % Borrowing constraint
    Ny = size(Xw{1},2);  % Number of income states (perm & transitory)
    V = cell(1,N);       % Value function
    c = cell(1,N);       % Consumption policy
    % A. Retirement period (Nw+1 -> N)
    c{N} = Xr;
    V{N} = c{N}.^(1-gamma) / (1-gamma);
    for j = N-1:-1:Nw+1
        Va = (1+r) * c{j+1}.^(-gamma);
        cj = (beta*ps(j+1) * Va).^(-1/gamma);
        xj = cj + a;  % Endogenous cash-on-hand grid
        Vj = cj.^(1-gamma) / (1-gamma) + beta*ps(j+1) * V{j+1};
        if amin < xj(1)  % Add x = a(1), c = 0 point to interpolation set
            cj = interp1([amin;xj],[0;cj],Xr,method,'extrap');
        else
            cj = interp1(xj,cj,Xr,method,'extrap');
        end
        V{j} = interp1(xj,Vj,Xr,method,'extrap');
        c{j} = max(min(cj,Xr-amin),0);  % Enforce borrowing constraint
    end
    % B. Retirement age (Nw)
    Va = repmat((1+r) * c{Nw+1}.^(-gamma),1,Ny);
    cj = (beta*ps(Nw+1) * Va).^(-1/gamma);
    xj = cj + a;  % Endogenous cash-on-hand grid
    Vj = cj.^(1-gamma) / (1-gamma) + beta*ps(Nw+1) * repmat(V{Nw+1},1,Ny);
    for k = 1:size(cj,2)
        xjk = xj(:,k);
        cjk = cj(:,k);
        if amin < xj(1,k)  % Add x = a(1), c = 0 point to interpolation set
            xjk = [amin;xjk];
            cjk = [0;cjk];
        end
        cj(:,k) = interp1(xjk,cjk,Xw{Nw}(:,k),method,'extrap');
        Vj(:,k) = interp1(xj(:,k),Vj(:,k),Xw{Nw}(:,k),method,'extrap');
    end
    c{Nw} = max(min(cj,Xw{Nw}-amin),0);  % Enforce borrowing constraint
    V{Nw} = Vj;
    % C. Working-life period (1 -> Nw-1)
    for j = Nw-1:-1:1
        Va = (1+r) * c{j+1}.^(-gamma);
        cj = (beta*ps(j+1) * Va * P').^(-1/gamma);
        xj = cj + a;  % Endogenous cash-on-hand grid
        Vj = cj.^(1-gamma) / (1-gamma) + beta*ps(j+1) * V{j+1} * P';
        for k = 1:size(cj,2)
            xjk = xj(:,k);
            cjk = cj(:,k);
            if amin < xj(1,k)  % Add x = a(1), c = 0 point to interpolation set
                xjk = [amin;xjk];
                cjk = [0;cjk];
            end
            cj(:,k) = interp1(xjk,cjk,Xw{j}(:,k),method,'extrap');
            Vj(:,k) = interp1(xj(:,k),Vj(:,k),Xw{j}(:,k),method,'extrap');
        end
        c{j} = max(min(cj,Xw{j}-amin),0);  % Enforce borrowing constraint
        V{j} = Vj;
    end
end
