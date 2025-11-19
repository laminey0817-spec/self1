function estFin = EMBiGAMP(Y,problem,BiGAMPopt, EMopt)
% Y: the noisy data matrix 带噪矩阵Y
% estFin: Structure containing final BiG-AMP outputs 存的是计算结果

%% 统计循环次数
sum_bigamp = 0;

%% Handle Options
%Get problem dimensions
M = problem.M; % 天线数
L = problem.L; % 信号长度
N = problem.N; % 用户数

BiGAMPopt.tol = 1e-4; %largest allowed tolerance for a single EM iteration

%Ensure that EM outputs are calculated
BiGAMPopt.saveEM = 1;
%% Initial Setup

%Set initial noise variance
nuw = EMopt.noise_var; % 噪声方差

lambda = EMopt.lambda; % 稀疏性，由于此场景下是没有稀疏性的，故lambda为1，设置为1-1e-12是为了数值稳定

%Initialize loop
t = 0; % EM-BiGAMP的循环次数
stop = 0;

nuX = EMopt.nuX; % 信号方差（for 初始化）
nuA = EMopt.nuA; % 信道方差（for 初始化）

%Initialize xhat with zeros
xhat = zeros(N,L); % X均值的初始化

%Initialize Ahat with random
Ahat = (randn(M,N) + randn(M,N)*1i) / sqrt(2) .* sqrt(nuA); % 信道均值的初始化（随机）

%Set init
BiGAMPopt.xhat0 = xhat;             % 信号均值的初始值（全0）
BiGAMPopt.Ahat0 = Ahat;             % 信道均值的初始值（随机）
BiGAMPopt.Avar0 = ones(M,N) .* nuA; % 信号方差的初始值（全1）
BiGAMPopt.xvar0 = ones(N,L) .* nuX; % 信道方差的初始值（全1）

%% Main Loop
%EM iterations
while ~stop
    
    %Increment time exit loop if exceeds maximum time
    t = t + 1;
    
    %Prior on A
    gA = CAwgnEstimIn(0, EMopt.nuA);
    
    %Prior on X
    gXbase = QPSKEstim();
    gX = SparseEstim(gXbase,lambda); 

    %% Output log likelihood
    gOut = CAwgnEstimOut(Y, nuw);
        
    % Run BiG-AMP
    estFin2 = BiGAMP(gX, gA, gOut, BiGAMPopt);
    
    %Check for estimate tolerance threshold
    if t >= EMopt.maxEMiter
        stop = stop + 1;
    end
    
    %Reinitialize EM estimates
    % 重启：只保留上次EM循环的
    BiGAMPopt.Ahat0 = estFin2.Ahat;
    BiGAMPopt.Avar0 = estFin2.Avar;
    BiGAMPopt.step = BiGAMPopt.stepMin; 
    
    % 统计bigamp循环次数
    sum_bigamp = sum_bigamp + estFin2.BiGcycletimes;

end

%% Cleanup

%Update finals
estFin = estFin2;

% 记录循环次数
estFin.EMcycletimes = t;
estFin.BiGcycletimes = sum_bigamp;
