function [out, share_bank, share_bond, share_none] = calc_Bi_BLP(Par, inp, f, Er, flag)

% This fucntion calcualtes Bi given ri when flag = 'r'
%            It calculates the inverse demand ri given Bi when flag = 'B'

% f is the federal funds rate
% A is the demand shifter
% Er is the estimated loan rate charged by other banks

% [inp, f, A, Er] should contain vestors of the same size and/or scalers

    max_length =  max([length(f), length(Er), length(inp)]);
    if (length(f)  == max_length | length(f)   ==1)...
      &(length(Er)  == max_length | length(Er)  ==1)...
      &(length(inp) == max_length | length(inp) ==1);
    else
        disp('matrix dimensions do not match in calc_Bi_BLP')
    end
       
    Ef = interp1(Par.fGrid, Par.Ef_Grid, f);
    
    exp_corporate  = exp(Par.nuC - Par.alpha * (Ef + 0*Par.meanA + (Ef-Par.meanf)*Par.rhoAf1) + 0);
    exp_0          = 1;
    exp_Ebank      = exp(Par.nuB - Par.alpha * Er)*(Par.NB-1);

    if flag == 'r' 
        ri = inp;
        exp_banki = exp(Par.nuB - Par.alpha * ri);
        Bi = exp_banki./(exp_banki + exp_Ebank + exp_corporate + exp_0);
        out = Bi;
    end

    if flag == 'B' 
        Bi = inp;
        exp_banki = (exp_Ebank + exp_corporate + exp_0).*Bi./(1-Bi);
        ri = (Par.nuB - log(exp_banki))/Par.alpha;
        ri(Bi >= 1)  = nan;
        ri(Bi <= 0) = nan;
        out  = ri; 
    end

    share_bank = (exp_banki + exp_Ebank)./(exp_banki + exp_Ebank + exp_corporate + exp_0);
    share_bond = exp_corporate./(exp_banki + exp_Ebank + exp_corporate + exp_0);
    share_none = exp_0./(exp_banki + exp_Ebank + exp_corporate + exp_0);

end