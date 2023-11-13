addpath('Script')
addpath('Library')
addpath('Data')

for iter = [1,2,3]
    if iter == 1
        disp(['============== Full Sample ==============='])
        st_dataset=readtable('Data\net_charge_off_raw_1.csv');
    elseif iter ==2 
        disp(['============== Early Sample ==============='])
        st_dataset=readtable('Data\net_charge_off_raw_2.csv');
    elseif iter ==3
         disp(['============== Late Sample ==============='])
         st_dataset=readtable('Data\net_charge_off_raw_3.csv');
    end
    columnname = st_dataset.Properties.VariableNames;

    for i = 1:length(columnname)
        eval([char(columnname(i)),'= table2array(st_dataset(:,i));']);
    end

    id         = rssdid;
    ffr        = FEDFUNDS;
    chargeoff1 = netloanchargeoffs_loans;
    chargeoff0 = netloanchargeoffs_loans_L1y;
    share      = assets_share;
    year_old = year;
    
    [unique_year,unique_yindex] = unique(year);
    [unique_id,unique_idindex] = unique(id);
    unique_ffr     = ffr(unique_yindex); 


    %% Calibrate the process of FFR

    cutoff0f = sum(prctile(ffr,[1:100]) < 0.001*100); % The lower bound in our simulation is 0.001
    cutoff1f = 100-sum(prctile(ffr,[1:100]) > 0.08*100); % The upper bound in our simulation is 0.08

    ffr = winsor(ffr,[cutoff0f cutoff1f]);
    muf = mean(log(ffr));
    stdf = std(log(ffr));

    std_muf = stdf/sqrt(length(unique_year));
    std_stdf = sqrt(2*stdf^4/(length(unique_year)-1));

    ffr0 = ffr*0;
    for i = 1:length(unique_year)
        ffr0 = ffr0 + unique_ffr(i)*(year == unique_year(i)+1);
    end

    dummy = find((share >=0).*(year>unique_year(1)));
    ffr        = ffr(dummy);
    ffr0      = ffr0(dummy);
    chargeoff1 = chargeoff1(dummy);
    chargeoff0 = chargeoff0(dummy);
    share      = share(dummy);
    year   = year(dummy);
    id = id(dummy);

    y = log(ffr)-muf;
    x = log(ffr0)-muf;

    [betaf]= lscov(y,x);
    std_betaf =  (mean(x.*x))^(-1) * mean (x.*mean((y-betaf*x).^2).*x) * (mean(x.*x))^(-1)/sqrt(length(unique_year));

    str   = ['mu_ffr= ', num2str(muf)]; disp(str)
    str   = ['std_ffr= ', num2str(stdf)]; disp(str)
    str   = ['rho_ffr= ', num2str(betaf)]; disp(str)




    %% Calibrate Chargeoff process

    cutoff0A = 1;
    cutoff1A = 99;

    chargeoff0 = winsor(chargeoff0,[cutoff0A cutoff1A]);

    dummy = (~isnan(chargeoff0)).*(~isnan(chargeoff1))>0;
    chargeoff1 = chargeoff1(dummy);
    chargeoff0 = chargeoff0(dummy);
    ffr = ffr(dummy);
    ffr0 = ffr0(dummy);
    share = share(dummy);
    id = id(dummy);

    y = log(chargeoff1(chargeoff1>0));
    x = y*0+1;
    w = winsor(share(chargeoff1>0),[5 95]);

    muA= lscov(x,y,w);
    std_muA =  (mean(x.*w.*w.*x))^(-1) * mean (x.*w.*w*mean((y-muA*x).^2).*w.*w.*x) * (mean(x.*w.*w.*x))^(-1)/sqrt(length(unique_id));
    str   = ['mu_chargeoff= ',num2str(muA)]; disp(str)

    y = demean(id(chargeoff1>0), log(chargeoff1(chargeoff1>0)));
    y = winsor(y.^2,[1 99]);
    x = y*0+1;
    w = winsor(share(chargeoff1>0),[5 95]);

    stdA= sqrt(lscov(x,y,w));
    std_varA =  (mean(x.*w.*w.*x))^(-1) * mean (x.*w.*w*mean((y-stdA*x).^2).*w.*w.*x) * (mean(x.*w.*w.*x))^(-1)/sqrt(length(unique_id));
    std_stdA =  sqrt(std_varA )/2;
    str   = ['std_chargeoff= ',num2str(stdA)]; disp(str)



    %% Calibrate serial correlation and persistence

    dummy = find((chargeoff1 > 0) & (chargeoff0 > 0));
    ffr        = ffr(dummy);
    ffr0      = ffr0(dummy);
    chargeoff1 = chargeoff1(dummy);
    chargeoff0 = chargeoff0(dummy);
    share      = share(dummy);
    year   = year(dummy);
    id = id(dummy);

    w = winsor(share,[5 95]);
    y = log(chargeoff1); y = winsor(y - muA,[1 99]);
    x  = log(ffr)-muf;
    betafA= lscov(x,y,w);
    std_betafA =  (mean(x.*w.*w.*x))^(-1) * mean (x.*w.*w*mean((y-betafA*x).^2).*w.*w.*x) * (mean(x.*w.*w.*x))^(-1)/sqrt(length(unique_id));

    str   = ['corr_chargeoff_ffr= ', num2str(betafA)]; disp(str)


    %% The chargeoff properties
 
    ffr        = FEDFUNDS;
    chargeoff1 = netloanchargeoffs_loans;
    chargeoff0 = netloanchargeoffs_loans_L1y;
    share      = assets_share;
    year = year_old;

    ffr0 = ffr*0;
    for i = 1:length(unique_year)
        ffr0 = ffr0 + unique_ffr(i);
    end
   chargeoff0 = (chargeoff0 > 0).*chargeoff0 + chargeoff1.*(~(chargeoff0 > 0)).*(chargeoff1>0);

   dummy = find((chargeoff1 > 0).*(chargeoff0 > 0).*(share > 0).*(~isnan(chargeoff0)).*(~isnan(chargeoff1))>0);
    ffr        = ffr(dummy);
    ffr0      = ffr0(dummy);
    chargeoff1 = chargeoff1(dummy);
    chargeoff0 = chargeoff0(dummy);
    share      = share(dummy);

    log_chargeoff1 = log(chargeoff1) - (log(ffr)-muf)*betafA;
    log_chargeoff0 = log(chargeoff0) - (log(ffr)-muf)*betafA;

    w = winsor(share, [5,95]);
    y = winsor((log_chargeoff1) - muA, [1 99]);
    x = winsor((log_chargeoff0) - muA, [1 99]);

    betaA = lscov(x,y,w); 
    std_betaA =  (mean(x.*w.*w.*x))^(-1) * mean (x.*w.*w*mean((y-betaA*x).^2).*w.*w.*x) * (mean(x.*w.*w.*x))^(-1)/sqrt(length(unique_id));

    epsilonA = y - x*betaA; 
    str   = ['rho_chargeoff = ', num2str(betaA)]; disp(str)



    %% Save results
    temp = [muf, betaf, stdf, betafA, muA, betaA, stdA ;
        std_muf, std_betaf, std_stdf, std_betafA, std_muA, std_betaA, std_stdA];
    
    % Print parameter estimates
    if iter ==1
        fileID = fopen('Results\FFR_Chargeoff.txt', 'a');        
        fprintf(fileID, '\n\n FFR and chargeoff process \n');
        fprintf(fileID, '[muf, rhof, stdf, corrAf, muA, rhoA, stdA] \n\n');
        fprintf(fileID,'Full Sample\n');
     elseif iter ==2
        fileID = fopen('Results\FFR_Chargeoff.txt', 'a');
        fprintf(fileID,'Early Sample \n');
     elseif iter ==3
        fileID = fopen('Results\FFR_Chargeoff.txt', 'a');
        fprintf(fileID,'Late Sample \n');
    end

    fprintf(fileID,'%6.1f %6.1f %6.1f %6.2f %6.1f %6.1f %6.1f\n ', temp(1,:)');
    fclose(fileID);

    % Print standard error
    if iter ==1
        fileID = fopen('Results\standard error.txt', 'a');
        fprintf(fileID, '\n\n FFR and chargeoff process \n');
        fprintf(fileID, '[muf, rhof, stdf, corrAf, muA, rhoA, stdA] \n\n');
        fprintf(fileID,'Full Sample\n');
    elseif iter ==2
        fileID = fopen('Results\standard error.txt', 'a');
        fprintf(fileID,'Early Sample \n');
    elseif iter ==3
        fileID = fopen('Results\standard error.txt', 'a');
        fprintf(fileID,'Late Sample \n');
    end
    fprintf(fileID,'%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n ', temp(2,:)');
    fclose(fileID);
    
end
 
