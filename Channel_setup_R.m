clc;clear;close all;

%% 参数设置
Tx = 4; Rx = 4;  % 天线配置4T4R
N_SC = 48;  % 子载波数(4RB)
N_Sym = 14;  % 每帧OFDM符号数
SeedRange = 50; % 生成的信道总个数= 13*SeedRange


carrier = nrCarrierConfig;  % 用户信号子载波配置
carrier.SubcarrierSpacing = 15;
carrier.NSizeGrid = 52;

grid_Tx = nrResourceGrid(carrier,Tx);  % 资源网格维度
grid_Tx = randn(size(grid_Tx));

[waveform_T,info_OFDM] = nrOFDMModulate(carrier,grid_Tx);  % 用户信号OFDM调制得到的时域波形

%% TDL信道
% 矩阵第1行表示TDLA30，第2行表示TDLB100，第3行表示TDLC300
PathDelays = 1e-9.*[0, 10, 15, 20, 25, 50, 65, 75, 105, 135, 150, 290
    0, 10, 20, 30, 35, 45, 55, 120, 170, 245, 330, 480
    0, 65, 70, 190, 195, 200, 240, 325, 520, 1045, 1510, 2595];
AveragePathGains = [-15.5, 0, -5.1, -5.1, -9.6, -8.2, -13.1, -11.5, -11, -16.2, -16.6, -26.2
    0, -2.2, -0.6, -0.6, -0.3, -1.2, -5.9, -2.2, -0.8, -6.3, -7.5, -7.1
    -6.9, 0, -7.7, -2.5, -2.4, -9.9, -8, -6.6, -7.1, -13, -14.2, -16];
% DelaySpread = 1e-9.*[30 100 300];
DopplerShift = [5,400,100];

% 选择信道TDLA30(k=1),TDLB100(k=2),TDLC300(k=3)
k = 1;

R_1_sum = zeros(48,48,Rx);
R_2_sum = zeros(48,48,Rx);
R_3_sum = zeros(48,48,Rx);
R_4_sum = zeros(48,48,Rx);

for slot = 1:SeedRange
    slot
    
    % 物理信道tdl
    tdl = nrTDLChannel;
    tdl.SampleRate = 30.72e6;
    tdl.DelayProfile = 'custom';
    %     tdl.DelaySpread = DelaySpread(k);  % 属性与 System object 的此配置无关
    tdl.PathDelays = PathDelays(k,:);
    tdl.AveragePathGains = AveragePathGains(k,:);
    tdl.MaximumDopplerShift = DopplerShift(k);  % 经典模型多普勒频偏取值TDLA30-5/10,TDLB100-400,TDLC300-100
    tdl.NumTransmitAntennas = Tx;
    tdl.NumReceiveAntennas = Rx;
    tdl.Seed = randi([1,SeedRange]);
    % tdl.Seed = 74;
    
    % 获取用户信道信息
    Info_TDLchannal0 = info(tdl);
    pathFilters = getPathFilters(tdl);
    [waveform_R,pathGains,sampleTimes] = tdl(waveform_T);  %信号经过信道tdl衰落
    
    hest = nrPerfectChannelEstimate(carrier,pathGains,pathFilters);  % 频域信道矩阵hest：N_SC-by-N_SYM-by-N_R-by-N_T
    
    for ii = 1:13
        H_part = hest(48*(ii-1)+1:48*ii,:,:,:); % 48 14 4 4
        H(:,:,:,:,13*(slot-1)+ii) = H_part;  % 生成的信道矩阵保存为H
        
        H11 = H_part; % 48 14 4 4

        % 用户1
        H_user1 = squeeze(H11(:,:,:,1)); % 48 14 4
        
        % 用户2
        H_user2 = squeeze(H11(:,:,:,2)); % 48 14 4
        
        % 用户3
        H_user3 = squeeze(H11(:,:,:,3));
        
        % 用户4
        H_user4 = squeeze(H11(:,:,:,4));
        
        % 计算频域自相关阵
        % 用户1的自相关矩阵
        R_1_temp = pagemtimes(H_user1,'none',H_user1,'ctranspose')/14; % 48 48 4
        % 不转置 * 共轭转置
        % 逐天线逐用户
        
        % 用户2的自相关矩阵
        R_2_temp = pagemtimes(H_user2,'none',H_user2,'ctranspose')/14;
        
        % 用户3的自相关矩阵
        R_3_temp = pagemtimes(H_user3,'none',H_user3,'ctranspose')/14;
        
        % 用户4的自相关矩阵
        R_4_temp = pagemtimes(H_user4,'none',H_user4,'ctranspose')/14;
        
        % 求平均
        R_1_sum = R_1_sum + R_1_temp;
        R_2_sum = R_2_sum + R_2_temp;
        R_3_sum = R_3_sum + R_3_temp;
        R_4_sum = R_4_sum + R_4_temp;
    end
end

R_1 = R_1_sum/(13 * SeedRange); % 48 48 4
R_2 = R_2_sum/(13 * SeedRange);
R_3 = R_3_sum/(13 * SeedRange);
R_4 = R_4_sum/(13 * SeedRange);


%% 统计TDLA B C的方差
RB = 4;
J = 12;
variance_H = zeros(size(H,5),1);
for idx = 1:size(H,5)
    
    H1 = H(:,:,:,:,idx); % 48 14 4 4
    
    %     i_RB = randperm(RB,1);
    variance_H_temp = zeros(1,RB);
    for i_RB = 1:RB
        H11 = H1( J * (i_RB-1) + 1 : J * i_RB ,[4:11,13 14],:,:);
        % 第i_RB个RB，除去导频位置 维度：12 10 4 4
        H111 = permute(H11,[3,4,1,2]); % 4 4 12 10
        
        deta_H = (H111(:,:,J,1) - H111(:,:,1,1))/J;
        
        % input variance of deta_H： sum(abs(deta_H(:)).^2) / numel(deta_H)
        variance_H_temp(i_RB) =  sum(abs(deta_H(:)).^2) / numel(deta_H);
    end
    variance_H(idx) = sum(variance_H_temp)/RB;
end

average_variance_H = sum(variance_H)/numel(variance_H);

if k == 1
    average_variance_H_A = average_variance_H;
    save('Channel_TDLA30_R.mat','H','R_1','R_2','R_3','R_4','average_variance_H_A');
elseif k == 2
    average_variance_H_B = average_variance_H;
    save('Channel_TDLB100_R.mat','H','R_1','R_2','R_3','R_4','average_variance_H_B');
else
    average_variance_H_C = average_variance_H;
    save('Channel_TDLC300_R.mat','H','R_1','R_2','R_3','R_4','average_variance_H_C');
end
