function [cnt,list_X,list_Xvar,list_Ahat,list_Avar] = checksame(optIn,temp_X,temp_Xvar,temp_Ahat,...
    temp_Avar,cnt,list_X,list_Xvar,list_Ahat,list_Avar)

Priority_proportion = optIn.Priority_proportion;           
l = size(list_X,1); %list里的包数
cnt_new = cnt;
for li = 1:l
    if sum(temp_X == list_X(li,:)) >= Priority_proportion * length(temp_X)
        cnt_new(li) = cnt_new(li)+1; % X对应索引+1
        
        % 信道高斯合并
        Ahat1 = list_Ahat(:,li);
        Ahat2 = temp_Ahat;
        Avar1 = list_Avar(:,li);
        Avar2 = temp_Avar;
        gain = Avar1./(Avar1+Avar2);
        list_Ahat(:,li) = gain.*(Ahat2-Ahat1)+Ahat1;
        list_Avar(:,li) = gain.*Avar2;
        
        % 如果新的信号的方差<list中对应的，则替换掉list中的，否则保留
        if sum(temp_Xvar,'all') < sum(list_Xvar(li,:))
            list_X(li,:) = temp_X; 
            list_Xvar(li,:) = temp_Xvar;  
        end
 
        break;
    end
end

if isequal(cnt_new,cnt)  % 此包没有和list中的任何一个包相同，就添加到list中
    list_X = [list_X;temp_X];
    list_Xvar = [list_Xvar;temp_Xvar];    
    cnt_new = [cnt_new;1];
    list_Ahat = [list_Ahat,temp_Ahat];
    list_Avar = [list_Avar,temp_Avar];
end

cnt = cnt_new;

end


