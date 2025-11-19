function results = eliminate_permutation(optIn,estFinTemp)
    L = optIn.L;
    alphabet = optIn.alphabet;
    Xhat1 = estFinTemp.xhat;
    Ahat = estFinTemp.Ahat;
    Avar = estFinTemp.Avar;
    xvar = estFinTemp.xvar;      
    Nhat = size(Xhat1,1);
    
    % 消除相位模糊
    for i = 1:Nhat
        %First symbol is alphabet(1)
        [~,alphabet_loc] = min(abs(Xhat1(i,1) ./ alphabet - alphabet(1)));
        scale = alphabet(alphabet_loc);
        Xhat1(i,:) = Xhat1(i,:) ./ scale;
        Ahat(:,i) = Ahat(:,i) .* scale;
    end
    
    Xhat2 = Xhat1;
    % 硬判决
    for i = 1:Nhat
        for j = 1:L
            if abs(Xhat2(i,j)) > 0
                [~,loc] = min(abs(alphabet - Xhat2(i,j)));
                Xhat2(i,j) = alphabet(loc);
            end
        end
    end
     
    %% 消除相位模糊,硬判决后的数据
    results.xhat = Xhat2;
    results.Ahat = Ahat;
    results.xvar = xvar;
    results.Avar = Avar;
    
end