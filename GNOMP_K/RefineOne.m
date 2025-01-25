function [omega_hat,g_hat,Tau4_r] = RefineOne(Tau4,sigma,omega_hat,g_hat,Phy_matrix,R_s)
% 对单个目标进行牛顿修正

% 输入：
%   Tau4：等效阈值，三个维度的大小分别为：N(M)，T，4
%         其中第三个维度顺序为（实部上门限，实部下门限，虚部上门限，虚部下门限）
%   sigma：噪声方差，是一个标量
%   omega_hat：频率估计值，是一个标量
%   g_hat：幅度估计值，是一个标量
%   Phy_matrix：压缩矩阵，大小为[M,N]，非压缩场景则为单位矩阵
%   R_s：NOMP步骤单循环次数

% 输出：
%   omega_hat：修正后的频率估计值，是一个标量
%   g_hat：修正后的幅度估计值，是一个标量
%   Tau4_r：更新的等效阈值
for i = 1:R_s
    omega_hat = ref_freq(Tau4, omega_hat, g_hat, sigma, Phy_matrix);
    [g_hat,~, Tau4_r] = ref_amp(Tau4, omega_hat, g_hat, sigma, Phy_matrix);
end
end

