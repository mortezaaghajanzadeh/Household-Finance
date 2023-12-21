function X = makeitpoly(X,polydeg)

%%% this is just a simple handler to make polynomials terms of arrays
%%% inside a matrix


nParams  = size(X,2) ;
      
        

%% Terms of 2nd degree
    iCol = nParams ;
    if polydeg>=2

            for i=1:nParams
    
    
                % vectorized code
                nCol = numel(i:nParams) ;
    
                ColIdx = iCol+1:iCol+nCol ;
    
                X(:,ColIdx) = X(:,i).*X(:,i:nParams) ;
    
                iCol = iCol+nCol ;
    
        %         % non-vectorized code
        %         for ii=i:nParams
        % 
        %             iCol = iCol + 1 ;
        % 
        %             X(:,iCol) = X(:,i).*X(:,ii) ;
        % 
        %         end
    
            end


        if polydeg>=3
%% Terms of 3rd degree

            x0 = nParams+1;
            x2 = iCol ;
    
            for i=1:nParams
    
                % vectorized code
                nCol = numel(x0:x2) ;
    
                ColIdx = iCol+1:iCol+nCol ;
    
                X(:,ColIdx) = X(:,i).*X(:,x0:x2) ;
    
                iCol = iCol+nCol ;
    
                x0 = x0+nParams-(i-1) ; 
    
        %         % non-vectorized code        
        %         for ii=x0:x2
        % 
        %             iCol = iCol + 1 ;
        % 
        %             X(:,iCol) = X(:,i).*X(:,ii) ;
        % 
        %         end
        % 
        %         x0 = x0+nParams-(i-1) ;     
    
            end

            if polydeg>=4
                error('this handler is not desgin for powers of more than 3!')
            end

        end

    end



%% Intercept
if polydeg>=1

    iCol = iCol +1 ;
    X(:,iCol) = 1 ;

end



end

