function [OmegaList,XList] = RefineAll(OmegaList,XList,Tau4_r,Phy_matrix,sigma,R_s,R_c)
% 对所有已经估计出的目标进行循环牛顿修正

% 输入：
%   OmegaList：频率列向量
%   XList：幅度列向量（元素为复数）
%   Tau4_r：等效阈值，三个维度的大小分别为：N(M)，T，4
%         其中第三个维度顺序为（实部上门限，实部下门限，虚部上门限，虚部下门限）
%   Phy_matrix：压缩矩阵，大小为[M,N]，非压缩场景则为单位矩阵
%   sigma：噪声方差，是一个标量
%   R_s：NOMP步骤单循环次数
%   R_c：NOMP循环修正步骤的循环次数

% 输出：
%   OmegaList：修正后的频率列向量
%   XList：修正后的幅度列向量（元素为复数）

N = size(Phy_matrix,2);
Khat = length(OmegaList);
for idx_R = 1:R_c
    order = 1:Khat;
    for j_tar = 1:Khat
        l = order(j_tar);
        omega_l = OmegaList(l);
        gain_l = XList(l,:);
        z_l = Phy_matrix*exp(1j*(0:N-1)'*omega_l)*gain_l;
        y_l = Tau4_r+cat(3,real(z_l),real(z_l),imag(z_l),imag(z_l));
        [omega_l,gain_l,Tau4_r] = RefineOne(y_l,sigma,omega_l,gain_l,Phy_matrix,R_s);
        OmegaList(l) = omega_l;
        XList(l,:) = gain_l;
    end
end
end

