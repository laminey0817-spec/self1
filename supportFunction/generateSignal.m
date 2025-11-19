function [Signal,optIn] = generateSignal(optIn, H,idx)
M = optIn.M;
N = optIn.N;
L = optIn.L;
RB = optIn.RB;
alphabet = optIn.alphabet;
%% 产生信号
for ii = 1:optIn.RB
    
    X_bits = randi([0,3],optIn.N, optIn.L);
    X_R1 = pskmod(X_bits,4,0,'gray'); % 格雷码
 
    % 加上一个参考符号
    X_R1(:,1) = alphabet(1); 
    X_bits(:,1) = 0;
        
    % 网格的处理（符号）
    X_R11 = reshape(X_R1,[optIn.N,optIn.L/optIn.length_t,optIn.length_t]); % 4 6 10 变到网格里
    optIn.X_allRB(:,(ii-1)*optIn.L/optIn.length_t+1:ii*optIn.L/optIn.length_t,1:optIn.length_t) = X_R11; 
    % 4 48 10（变回整个大网格中）
    optIn.X_allRB1(:,:,ii) = X_R1; % 4 60 8 或者 4 120 4 
    
    % 对比特做转换
    X_bits_grid = reshape(X_bits,[optIn.N,optIn.L/optIn.length_t,optIn.length_t]); % 4 6 10 变到网格里
    optIn.X_bits_allRB(:,(ii-1)*optIn.L/optIn.length_t+1:ii*optIn.L/optIn.length_t,1:optIn.length_t) = X_bits_grid;
    % 4 48 10（变回整个大网格中）
    
end

%% 信道 + 接收信号
X_allRB1 = optIn.X_allRB1;
H1 = H(:,:,:,:,idx); % 48 14 4 4 取出所有信道的第几个
length_f = optIn.L/optIn.length_t; % 每次bigamp的实现所占子载波个数
Y_nonoise = zeros(M,L,RB);
H11_g = zeros(M,N,L,RB);
for ii = 1:RB
    
    H11 = H1((ii-1)*length_f+1:ii*length_f,[4:11,13 14],:,:); % 6 10 4 4
    H11_p = permute(H11,[3,4,1,2]); % 4 4 6 10
    H_real(:,:,(ii-1)*length_f+1:ii*length_f,1:10) = H11_p; % 4 4 6 10
    H11_g(:,:,:,ii) = reshape(H11_p,[M,N,L]); % 4 4 60 8 
    
    for iL = 1:L
        Y_nonoise(:,iL,ii) = H11_g(:,:,iL,ii) * X_allRB1 (:,iL,ii); % 接收信号Y
    end
    
end

%% 函数输出变量
Signal.Y_nonoise = Y_nonoise; % 4 60 8 或者 4 120 4
Signal.H_real = H_real; % 4 4 48 10
end