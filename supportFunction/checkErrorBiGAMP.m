function  errRes = checkErrorBiGAMP(optIn,result_all1,results2,Signal)

N = optIn.N;
M = optIn.M;
alphabet = optIn.alphabet; % 星座点
X_allRB = optIn.X_allRB; % 4 48 10 真实

X_wiener_zhuzaibo = results2.X_wiener_zhuzaibo; % 维纳后的 没有硬判决
X_junzhixielv_wiener = results2.X_junzhixielv_wiener;

H_wiener_white = results2.H_wiener_white; % 4 4 48
H_wiener_nowhite = results2.H_wiener_nowhite; % 4 4 48

RB = optIn.RB;
length_f = optIn.L/10; % 12
length_t = optIn.length_t; % 10

%% MSE（逐载波修正）
H_real = Signal.H_real; % 4 4 48 10 真实
Ahat_final = result_all1.Ahat; % 4 4 8 估计的

% 第一次big
H_1 = zeros(4,4,48,10);
for ii = 1:RB
    H_1( :,:,(ii-1)*length_f+1:ii*length_f,1:10) = Ahat_final(:,:,ii).* ones(4,4,length_f,length_t);
end
MSE1 = sum(abs(H_real - H_1).^2,'all');

% 白化+维纳滤波
for ij = 1:48
    for ii = 1:length_t
        H_2(:,:,ij,ii) = H_wiener_white(:,:,ij);
    end
end
MSE2 = sum(abs(H_real - H_2).^2,'all');

% 不白化+维纳滤波
for ij = 1:48
    for ii = 1:length_t
        H_3(:,:,ij,ii) = H_wiener_nowhite(:,:,ij);
    end
end
MSE3 = sum(abs(H_real - H_3).^2,'all');

%% SER
% 第一次的结果
Xhat_grid_eliminat = result_all1.Xhat_grid_eliminat; % 4 48 10

SER1 = zeros(N,1);
for i = 1:N
    a = squeeze(Xhat_grid_eliminat(i,:,:) - X_allRB(i,:,:));
    SER1(i) = sum(abs(a)>1e-3,'all');
end

% 第二次的结果（逐载波修正）
% 硬判决
Xhat2_hard = X_wiener_zhuzaibo; % 4 48 10
for i1 = 1:N
    for i2 = 1:48
        for i3 = 1:10
            [~,loc] = min(abs(alphabet - Xhat2_hard(i1,i2,i3)));
            Xhat2_hard(i1,i2,i3) = alphabet(loc);
        end
    end
end

% 修正后的误码率（逐载波修正）
SER2 = zeros(N,1);
for i = 1:N
    a = squeeze(Xhat2_hard(i,:,:) - X_allRB(i,:,:));
    SER2(i) = sum(abs(a)>1e-3,'all');
end

% 均值斜率
% 硬判决
Xhat2_hard_junzhixielv = X_junzhixielv_wiener; % 4 48 10
for i1 = 1:N
    for i2 = 1:48
        for i3 = 1:10
            [~,loc] = min(abs(alphabet - Xhat2_hard_junzhixielv(i1,i2,i3)));
            Xhat2_hard_junzhixielv(i1,i2,i3) = alphabet(loc);
        end
    end
end

% 修正后的误码率（逐载波修正）
SER2_junzhixielv = zeros(N,1);
for i = 1:N
    a = squeeze(Xhat2_hard_junzhixielv(i,:,:) - X_allRB(i,:,:));
    SER2_junzhixielv(i) = sum(abs(a)>1e-3,'all');
end


%% BER

% 第一次的结果
% 解调Xhat_grid_eliminat回比特流
for ii = 1:10
    Xhat_bits1(:,:,ii) = pskdemod(Xhat_grid_eliminat(:,:,ii),4,0,'gray'); % big
end

% 换成二进制
error_bit1 = 0;
for Ni = 1:N
    for fi = 1:48
        for li = 1:10
        a = dec2bin(Xhat_bits1(Ni,fi,li),2);
        b = dec2bin(optIn.X_bits_allRB(Ni,fi,li),2);
        c = a - b;
        temp = sum(abs(c),'all');
        error_bit1 = error_bit1 + temp;
        end
    end
end

% 第二次的结果（逐载波修正）
% 解调Xhat2_hard回比特流
for ii = 1:10
    Xhat_bits2(:,:,ii) = pskdemod(Xhat2_hard(:,:,ii),4,0,'gray'); % big
end

% 换成二进制
error_bit2 = 0;
for Ni = 1:N
    for fi = 1:48
        for li = 1:10
        a = dec2bin(Xhat_bits2(Ni,fi,li),2);
        b = dec2bin(optIn.X_bits_allRB(Ni,fi,li),2);
        c = a - b;
        temp = sum(abs(c),'all');
        error_bit2 = error_bit2 + temp;
        end
    end
end


% 第二次的结果（均值斜率）
% 解调Xhat2_hard回比特流
for ii = 1:10
    Xhat_bits2_junzhixielv(:,:,ii) = pskdemod(Xhat2_hard_junzhixielv(:,:,ii),4,0,'gray'); % big
end

% 换成二进制
error_bit2_junzhixielv = 0;
for Ni = 1:N
    for fi = 1:48
        for li = 1:10
        a = dec2bin(Xhat_bits2_junzhixielv(Ni,fi,li),2);
        b = dec2bin(optIn.X_bits_allRB(Ni,fi,li),2);
        c = a - b;
        temp = sum(abs(c),'all');
        error_bit2_junzhixielv = error_bit2_junzhixielv + temp;
        end
    end
end

%%
errRes.SER1 = sum(SER1,'all') ./ (N*480); % 第一次big的误码率
errRes.SER2 = sum(SER2,'all') ./ (N*480); % 修正后的误码率
errRes.SER2_junzhixielv = sum(SER2_junzhixielv,'all') ./ (N*480); % 修正后的误码率

errRes.BER1 = error_bit1 ./ (N*480*2); % 修正后的误比特率
errRes.BER2 = error_bit2 ./ (N*480*2); % 修正后的误比特率
errRes.BER2_junzhixielv = error_bit2_junzhixielv ./ (N*480*2); % 修正后的误比特率

errRes.MSE1 = MSE1 ./ (N*M*480); % 第一次big的MSE
errRes.MSE2 = MSE2 ./ (N*M*480); % 白化+维纳滤波的MSE
errRes.MSE3 = MSE3 ./ (N*M*480); % 不白化+维纳滤波的MSE

end
