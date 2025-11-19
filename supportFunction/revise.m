function  [results2,result_all1] = revise(optIn,result_all1)
N = optIn.N;
M = optIn.M;
L = optIn.L;
J = optIn.J;
RB = optIn.RB; % 把48*10的信道分成了多少块
Y_allRB = optIn.Y_allRB; % 4 120 4 或者 4 60 8
X_allRB = optIn.X_allRB1; % 4 120 4 或者 4 60 8
threshold_var_revise = optIn.threshold_var_revise; % 认为大于此方差则是有问题的解，则进行修正
sigma2 = optIn.Wvar; 
X_ML_matrix = optIn.X_ML_matrix; % for ML，X的所有中可能 

length_f = optIn.L/10; % 每次实现对应的子载波数
length_t = optIn.length_t; % 每次实现对应的OFDM symbols
xhat = result_all1.xhat; % 第一次big的结果（均值）4 120 4
xvar = result_all1.xvar; % 第一次big的结果（方差）

% 自相关矩阵
R_1 = optIn.R_1;
R_2 = optIn.R_2;
R_3 = optIn.R_3;
R_4 = optIn.R_4;

%% 找到排列的顺序
compareX_Xhat = zeros(N,N,RB);
for ii = 1:RB
    for i = 1:N
        for i1 = 1:N
            compareX_Xhat(i,i1,ii) =  sum(X_allRB(i,:,ii) ~= xhat(i1,:,ii));
            % 第i行第j列表示：X的第i个包和Xhat的第j个包相差的个数
        end
    end
end

for ij = 1:RB
    compareX_Xhat_temp = compareX_Xhat(:,:,ij);
    for i = 1:N
        [n,~] = min(compareX_Xhat_temp,[],2);
        [~,p1] = min(n);
        [~,p2] = min(compareX_Xhat_temp(p1,:));
        
        P(i,:) = [p1,p2];
        compareX_Xhat_temp(p1,:) = L;
        compareX_Xhat_temp(:,p2) = L;
    end
    match_list = sortrows(P,1);
    positionH_X(:,ij)  = match_list(:,2); % 用户排列模糊性的矩阵
end

%% 消除排列的模糊性
% 消除排列模糊性后的信号矩阵
for ii = 1:RB
    for ni = 1:N
        xhat_new(ni,:,ii) = xhat(positionH_X(ni,ii),:,ii); % 4 120 4 或者 4 60 8
        xvar_new(ni,:,ii) = xvar(positionH_X(ni,ii),:,ii); %（120和60是信号长度，4和8是把48*10的4RB分成了几块来分别实现）
    end
end

%% 拼回网格（消除排列的模糊性）

% 把每块拼回网格中
for ii = 1:RB
    xhat_grid(:,:,:,ii) =  reshape(xhat_new(:,:,ii),[N,length_f,length_t]);  % 4 12 10 4
    xvar_grid(:,:,:,ii) =  reshape(xvar_new(:,:,ii),[N,length_f,length_t]);  % 4 12 10 4
    Y_grid(:,:,:,ii) =  reshape(Y_allRB(:,:,ii),[N,length_f,length_t]); 
end

% 拼回48*10的4RB的大网格中
for ii = 1:RB
    Xhat_grid( :,(ii-1)*length_f+1:ii*length_f,1:10,:) = xhat_grid(:,:,:,ii); % 大网格
    Xvar_grid( :,(ii-1)*length_f+1:ii*length_f,1:10,:) = xvar_grid(:,:,:,ii); % 4 48 10
    Y_grid_big( :,(ii-1)*length_f+1:ii*length_f,1:10,:) = Y_grid(:,:,:,ii); % 4 48 10
end

%% 逐载波信道估计
for ii = 1:48 % 逐载波信道估计
    Y_sub = squeeze(Y_grid_big(:,ii,:)); % 4 10 逐载波的接收信号
    X_sub = squeeze(Xhat_grid(:,ii,:)); % 4 10 逐载波的Xhat
    Xvar_sub = squeeze(Xvar_grid(:,ii,:)); % 4 10 逐载波的Xvar
    
    %% 不白化，直接做信道估计
    [xvar_max,~] = max(Xvar_sub,[],1); % 不白化，故是对所有用户的xvar求和后找方差小的列
    [a,tag]=sort(xvar_max,'ascend');
    position_temp = tag(1:length_t*4/5);  % cnt升序排列的前8个对应的位置索引
    
    % LMMSE
    x1 = squeeze(X_sub(:,position_temp))';
    y1 = squeeze(Y_sub(:,position_temp))';
    Hhat_temp = x1' * pinv((x1 * x1') + sigma2 * eye(size(y1,1))) * y1;
    H_sub_nowhite(:,:,ii) = Hhat_temp'; % 4 4 48 得到逐载波的估计信道 一共48个载波
    
    %% 白化后再估计（未用此部分的结果）
    for iu = 1:N  
        % 找到对这个用户的可靠的部分
        Xvar_u = Xvar_sub(iu,:); % 10 白化，故是对单个用户的xvar求和后找方差小的列
        [a,tag]=sort(Xvar_u,'ascend');
        positon = tag(1:length_t*4/5);  % cnt升序排列的前8个对应的位置索引

        % 取出对于这个用户的可靠的部分
        X_sub_trust = X_sub(:,positon);
        Y_sub_trust = Y_sub(:,positon);
        
        % 白化
        X_u = X_sub_trust(iu,:).'; % 服务用户
        X_i = X_sub_trust;
        X_i(iu,:) = []; % 干扰用户
        X_I1 = X_i(1,:).'; % 其他三个干扰用户 % 8 1
        X_I2 = X_i(2,:).';
        X_I3 = X_i(3,:).';
        
        R_xx = 0.25 * (X_I1 * X_I1'+ X_I2 * X_I2'+ X_I3 * X_I3')+ sigma2 * eye(size(X_I1,1));
        P = pinv(sqrtm(R_xx));
        X_eq1 = (P * X_u).'; % 8 1
        Y_eq1 = (P * Y_sub_trust.').'; % 8 4
        
        % 白化后的等效噪声
        W = 1/(sum(abs(X_eq1).^2,'all')/length(X_eq1));
        % 逐用户的信道估计 LMMSE
        x1 = X_eq1';
        y1 = Y_eq1';
        Hhat_temp = x1' * pinv((x1 * x1') + W * eye(size(y1,1))) * y1;
        H_sub_wihte(:,iu,ii) = Hhat_temp'; % 4 4 48 逐用户的得到逐载波的估计信道
    end
end

%% 维纳滤波-无白化
for iu = 1:N
    % user1
    H_user1 = squeeze(H_sub_nowhite(:,1,:)).';  % 12 4
    W_p(:,:,iu) = squeeze(R_1(:,:,1))*pinv(squeeze(R_1(:,:,1))+(sigma2)*eye(size(R_1,1)));
    H_user1_Wiener(:,iu) = W_p(:,:,iu) * H_user1(:,iu);
    
    % user2
    H_user2 = squeeze(H_sub_nowhite(:,2,:)).';
    W_p(:,:,iu) = squeeze(R_2(:,:,2))*pinv(squeeze(R_2(:,:,2))+(sigma2)*eye(size(R_2,1)));
    H_user2_Wiener(:,iu) = W_p(:,:,iu) * H_user2(:,iu);
    
    % user3
    H_user3 = squeeze(H_sub_nowhite(:,3,:)).';
    W_p(:,:,iu) = squeeze(R_3(:,:,3))*pinv(squeeze(R_3(:,:,3))+(sigma2)*eye(size(R_3,1)));
    H_user3_Wiener(:,iu) = W_p(:,:,iu) * H_user3(:,iu);
    
    % user4
    H_user4 = squeeze(H_sub_nowhite(:,4,:)).';
    W_p(:,:,iu) = squeeze(R_4(:,:,4))*pinv(squeeze(R_4(:,:,4))+(sigma2)*eye(size(R_4,1)));
    H_user4_Wiener(:,iu) = W_p(:,:,iu) * H_user4(:,iu);
end

% 将逐天线的拼回
H_wiener_nowhite = zeros(M,N,48); % 4 4 48
H_wiener_nowhite(:,1,:) = H_user1_Wiener.';
H_wiener_nowhite(:,2,:) = H_user2_Wiener.';
H_wiener_nowhite(:,3,:) = H_user3_Wiener.';
H_wiener_nowhite(:,4,:) = H_user4_Wiener.';

%% 维纳滤波-有白化
for iu = 1:N
    % user1
    H_user1 = squeeze(H_sub_wihte(:,1,:)).';  % 12 4
    W_p(:,:,iu) = squeeze(R_1(:,:,1))*pinv(squeeze(R_1(:,:,1))+(sigma2)*eye(size(R_1,1)));
    H_user1_Wiener(:,iu) = W_p(:,:,iu) * H_user1(:,iu);
    
    % user2
    H_user2 = squeeze(H_sub_wihte(:,2,:)).';
    W_p(:,:,iu) = squeeze(R_2(:,:,2))*pinv(squeeze(R_2(:,:,2))+(sigma2)*eye(size(R_2,1)));
    H_user2_Wiener(:,iu) = W_p(:,:,iu) * H_user2(:,iu);
    
    % user3
    H_user3 = squeeze(H_sub_wihte(:,3,:)).';
    W_p(:,:,iu) = squeeze(R_3(:,:,3))*pinv(squeeze(R_3(:,:,3))+(sigma2)*eye(size(R_3,1)));
    H_user3_Wiener(:,iu) = W_p(:,:,iu) * H_user3(:,iu);
    
    % user4
    H_user4 = squeeze(H_sub_wihte(:,4,:)).';
    W_p(:,:,iu) = squeeze(R_4(:,:,4))*pinv(squeeze(R_4(:,:,4))+(sigma2)*eye(size(R_4,1)));
    H_user4_Wiener(:,iu) = W_p(:,:,iu) * H_user4(:,iu);
end

% 将逐天线的拼回
H_wiener_white = zeros(M,N,48); % 4 4 48
H_wiener_white(:,1,:) = H_user1_Wiener.';
H_wiener_white(:,2,:) = H_user2_Wiener.';
H_wiener_white(:,3,:) = H_user3_Wiener.';
H_wiener_white(:,4,:) = H_user4_Wiener.';

%% 逐载波的信号检测（修正）
X_wiener_zhuzaibo = Xhat_grid; % 4 48 10
for ii = 1:48
    
    % 取出某个子载波的信道和接收信号
    % 是否白化统计mse后基本一样，故目前使用不白化的方案
    H_sub = H_wiener_nowhite(:,:,ii); % 4 4 
    Y_sub = squeeze(Y_grid_big(:,ii,:)); % 4 10 逐载波的接收信号
    
    % 找到需要修正的列
    Xvar_sub = squeeze(Xvar_grid(:,ii,:)); % 4 10
    [xvar_max,~] = max(Xvar_sub,[],1);
    positon_X_revise = find(xvar_max > threshold_var_revise); % 这个子载波上需要修正的位置
    
    % ML
    if size(positon_X_revise,1) ~= 0 % 如果有方差>阈值的部分，则进行修正
        for ij = 1:length(positon_X_revise)
            % 取出这个子载波上需要修正的列
            Y_temp = Y_sub(:,positon_X_revise(ij)); % 4 1
            H_temp = H_sub; % 4 4
            
            % ML
            Y_HX = zeros(1,256);
            for im = 1:256
                Y_HX(im) = norm((Y_temp - H_temp * X_ML_matrix(:,im)),'fro');
            end
            
            [~,minY_HXlocation] = min(Y_HX);
            X_wiener_zhuzaibo(:,ii,positon_X_revise(ij)) = X_ML_matrix(:,minY_HXlocation);
            
        end
    end
end

% 均值斜率模型
X_wiener = junzhixielv(xhat_new,xvar_new,Y_grid_big,Xhat_grid,Xvar_grid,optIn);
results2.X_junzhixielv_wiener = X_wiener; % 均值斜率模型

% 逐载波
results2.X_wiener_zhuzaibo = X_wiener_zhuzaibo;
results2.H_wiener_white = H_wiener_white; % 4 4 48 10 有白化的
results2.H_wiener_nowhite = H_wiener_nowhite; % 4 4 48 10 无白化的

results2.positionH_X = positionH_X; % 排列模糊性的顺序
result_all1.Xhat_grid_eliminat = Xhat_grid;
result_all1.Xvar_grid_eliminat = Xvar_grid;
end


