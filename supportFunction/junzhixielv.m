function X_wiener = junzhixielv(Xhat,Xvar,Y_grid_big,Xhat_grid,Xvar_grid,optIn)

% Xhat Xvar 4 120 4 消除模糊性之后的
M = optIn.M;
N = optIn.N;
L = optIn.L;
J = optIn.J;
j0 = J/2;
Y_allRB = optIn.Y_allRB; % 4 120 4 或者 4 60 8
X_allRB = optIn.X_allRB;
RB = optIn.RB;
alphabet = optIn.alphabet;
length_t = optIn.length_t;
sigma2 = optIn.Wvar;
threshold_var_revise = optIn.threshold_var_revise;
X_ML_matrix = optIn.X_ML_matrix;

% 自相关矩阵
R_1 = optIn.R_1;
R_2 = optIn.R_2;
R_3 = optIn.R_3;
R_4 = optIn.R_4;

for ib = 1:RB
    
    Xhat_single = Xhat(:,:,ib);
    Xvar_single = Xvar(:,:,ib);
    Y = Y_allRB(:,:,ib);
    
    % 变回网格
    Xhat_single_grid = reshape(Xhat_single,[N,J,10]);  % 4 12 10
    Xvar_single_grid = reshape(Xvar_single,[N,J,10]);  % 4 12 10
    
    % 认为估计的非常准确的信号 -> 准确的估计信道
    [xvar_max,~] = max(Xvar_single,[],1);
    [a,tag]=sort(xvar_max,'ascend');
    positon_Xacc = tag(1:L/10);  % cnt升序排列的前12个对应的位置索引
    
    % 把X的均值和方差加上系数
    % 均值
    xhat_c_grid = zeros(M,J,length_t);
    for j = 1:J
        coefficient = j - J/2;
        xhat_c_grid(:,j,:) = coefficient * Xhat_single_grid(:,j,:);
    end
    xhat_c = reshape(xhat_c_grid,[M,L]);
    
    % 方差
    xvar_c_grid = zeros(M,J,length_t);
    for j = 1:J
        coefficient = (j - J/2)^2;
        xvar_c_grid(:,j,:) = coefficient * Xvar_single_grid(:,j,:);
    end
    xvar_c = reshape(xvar_c_grid,[M,L]);
    
    % 取对应索引位置的xhat和xvar，有系数的
    xhat_choose_c = xhat_c(:,positon_Xacc);
    xvar_choose_c = xvar_c(:,positon_Xacc);  %此处的xvar已经加上系数C了，不是全都是<threshold_var的
    
    % 取对应索引位置的xhat和xvar，无系数的
    xhat_choose = Xhat_single(:,positon_Xacc);

    % detaH 和 detaX的方差，转移到噪声方差上
    var_detaX = sum(xvar_choose_c,'all')/numel(xvar_choose_c);
    var_detaX_H = optIn.var_slope * var_detaX;

    Wnew = optIn.Wvar + var_detaX_H;
    Y_new0 = Y(:,positon_Xacc); % 4 12
    
    x_detax = [xhat_choose;xhat_choose_c]'; % 12 8
    Y_new = Y_new0'; %12 4
    
    % LMMSE
    cov = [ones(1,M),ones(1,M)*optIn.var_slope]; % 1 8
    cov_matrix = diag(cov,0); % 8 8
    delta_H_hermit = cov_matrix * x_detax' * ...
        pinv((x_detax * cov_matrix* x_detax') + Wnew * eye(size(Y_new,1))) * Y_new; % 8 4
    delta_H_hat = delta_H_hermit'; % 4 8
    
    % 信道均值
    H = delta_H_hat(:,1:M);
    % 信道斜率
    detaH = delta_H_hat(:,M+1:end);

    % 得到每个子载波上的信道
    H_f = zeros(M,N,J);
    for j1 = 1: J
        H_f(:,:,j1) = H + (j1- j0) * detaH; % 4 4 12
    end

    H_allsub(:,:,(ib-1)*J+1:ib*J) = H_f;
end

%% 维纳滤波-有白化
for iu = 1:N
    % user1
    H_user1 = squeeze(H_allsub(:,1,:)).';  % 12 4
    W_p(:,:,iu) = squeeze(R_1(:,:,1))*pinv(squeeze(R_1(:,:,1))+(sigma2)*eye(size(R_1,1)));
    H_user1_Wiener(:,iu) = W_p(:,:,iu) * H_user1(:,iu);
    
    % user2
    H_user2 = squeeze(H_allsub(:,2,:)).';
    W_p(:,:,iu) = squeeze(R_2(:,:,2))*pinv(squeeze(R_2(:,:,2))+(sigma2)*eye(size(R_2,1)));
    H_user2_Wiener(:,iu) = W_p(:,:,iu) * H_user2(:,iu);
    
    % user3
    H_user3 = squeeze(H_allsub(:,3,:)).';
    W_p(:,:,iu) = squeeze(R_3(:,:,3))*pinv(squeeze(R_3(:,:,3))+(sigma2)*eye(size(R_3,1)));
    H_user3_Wiener(:,iu) = W_p(:,:,iu) * H_user3(:,iu);
    
    % user4
    H_user4 = squeeze(H_allsub(:,4,:)).';
    W_p(:,:,iu) = squeeze(R_4(:,:,4))*pinv(squeeze(R_4(:,:,4))+(sigma2)*eye(size(R_4,1)));
    H_user4_Wiener(:,iu) = W_p(:,:,iu) * H_user4(:,iu);
end

% 将逐天线的拼回
H_wiener = zeros(M,N,48); % 4 4 48
H_wiener(:,1,:) = H_user1_Wiener.';
H_wiener(:,2,:) = H_user2_Wiener.';
H_wiener(:,3,:) = H_user3_Wiener.';
H_wiener(:,4,:) = H_user4_Wiener.';


%% 信号检测（修正）
X_wiener = Xhat_grid; % 4 48 10
for ii = 1:48
    
    % 取出某个子载波的信道和接收信号
    % 是否白化统计mse后基本一样，故目前使用不白化的方案
    H_sub = H_wiener(:,:,ii); % 4 4 
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
            X_wiener(:,ii,positon_X_revise(ij)) = X_ML_matrix(:,minY_HXlocation);
            
        end
    end
end


