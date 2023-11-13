function [demean_out,mean_out] = demean(ID, data)

    % This function demean each data series using fixed firm effect. 
    % out is the de-meaned version of data.
    
    if size(ID,1) == size(data,1)
    else disp('Sizes of the ID and data are inconssitent');
    end
    
    unique_ID = unique(ID);
    NID = length(unique_ID);
    Nob = length(ID);
    Ncol= size(data,2);
    
    mean_short= zeros(NID,Ncol);
    mean_long  = zeros(Nob,Ncol);
    
    for j = 1:NID
        indexj = find(ID==unique_ID(j));
        dataj   = data(indexj,:);
        
        nj = size(dataj,1);

        mean_short(j,:)= mean(dataj,1);
        mean_long(indexj,:)= repmat(mean_short(j,:),nj,1);

    end
    
    demean_out = data - mean_long;
    mean_out = mean_short;
    
end

