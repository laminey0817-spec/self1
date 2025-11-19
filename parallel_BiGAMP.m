function [results1,optIn] = parallel_BiGAMP(optIn)
setupPath;

maxTrials = optIn.maxTrials;
Y = optIn.Y;

%set BiGAMP option
opt = BiGAMPOpt();  % initialize the options object
opt.nit = optIn.inIt;  % limit iterations 内层循环次数
problem.M = optIn.M;
problem.L = optIn.L;
problem.N = optIn.N;
problem.spar = optIn.spar;

%set EMopt
EMopt = [];
EMopt.noise_var = optIn.Wvar;
EMopt.lambda = optIn.spar;
EMopt.maxEMiter = optIn.outIt; % BiG-AM外层循环次数
EMopt.alphabet = optIn.alphabet;
EMopt.nuX = 1; %初始化的时候的信号方差
EMopt.nuA = optIn.Avar; %初始化的时候的信道方差

%% BiGAMP
xhat_all = zeros(optIn.N,optIn.L,maxTrials);
xvar_all = zeros(optIn.N,optIn.L,maxTrials);
Ahat_all = zeros(optIn.M,optIn.N,maxTrials);
Avar_all = zeros(optIn.M,optIn.N,maxTrials);
Y_AX = zeros(1,maxTrials);
sum_EM_cycle = zeros(1,maxTrials);
sum_BiG_cycle = zeros(1,maxTrials);

for trial = 1:maxTrials
    
    %Run BiGAMP
    estFinTemp = EMBiGAMP(Y,problem,opt,EMopt);
    
    % 循环次数
    sum_EM_cycle(trial) = estFinTemp.EMcycletimes;
    sum_BiG_cycle(trial) = estFinTemp.BiGcycletimes;
        
    % 消除相位模糊，并做硬判决
    results = eliminate_permutation(optIn,estFinTemp);
    
    % 存进来的没有是 消除相位模糊,硬判决后的
    xhat_all(:,:,trial) = results.xhat; % 4 60 5
    xvar_all(:,:,trial) = results.xvar;
    Ahat_all(:,:,trial) = results.Ahat; % 4 4 5 
    Avar_all(:,:,trial) = results.Avar;
    
    % 计算Y-HX的值（消除相位模糊，硬判决后的，没有排序的）
    Y_AX(trial) = norm((Y - Ahat_all(:,:,trial) * xhat_all(:,:,trial)),'fro');
end

%% 优先队列
if optIn.Priority_queue  % 使用优先队列的方法
    cnt = ones(optIn.N,1);
    list_X = xhat_all(:,:,1);
    list_Xvar = xvar_all(:,:,1);
    list_Ahat = Ahat_all(:,:,1);
    list_Avar = Avar_all(:,:,1);
    
    for ti = 2:maxTrials
        for Ni = 1:optIn.N
            temp_X = xhat_all(Ni,:,ti);
            temp_Xvar = xvar_all(Ni,:,ti);
            temp_Ahat = Ahat_all(:,Ni,ti);
            temp_Avar = Avar_all(:,Ni,ti);
            [cnt,list_X,list_Xvar,list_Ahat,list_Avar] = checksame(optIn,temp_X,temp_Xvar,temp_Ahat,...
                temp_Avar,cnt,list_X,list_Xvar,list_Ahat,list_Avar);
        end
    end
    
    if size(cnt,1) > optIn.N*(maxTrials-1) || logical(sum(cnt>maxTrials))
        [~,minY_AXlocation] = min(Y_AX);
        results1.xhat_final = xhat_all(:,:,minY_AXlocation);
        results1.xvar_final = xvar_all(:,:,minY_AXlocation);
        results1.Ahat_final = Ahat_all(:,:,minY_AXlocation);
    else
        [rank,tag]=sort(cnt,'descend');
        positon = tag(1:optIn.N); % cnt降序排列的前N个对应的位置索引
        results1.xhat_final = list_X(positon,:);
        results1.xvar_final = list_Xvar(positon,:);
        results1.Ahat_final = list_Ahat(:,positon);
    end
    
else

    % 没有进行第二次BiG-AMP的结果（消除相位模糊,硬判决后的）
    [~,minY_AXlocation] = min(Y_AX);
    results1.xhat_final = xhat_all(:,:,minY_AXlocation);
    results1.xvar_final = xvar_all(:,:,minY_AXlocation);
    results1.Ahat_final = Ahat_all(:,:,minY_AXlocation);   
end

% 统计循环次数
results1.sum_EM_cycle = sum(sum_EM_cycle);
results1.sum_BiG_cycle = sum(sum_BiG_cycle);

end
