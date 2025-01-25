function [OmegaList,XList] = GNOMP_K(Tau4,sigma,Phy_matrix,zeta,tau,K,R_s,R_c,OverSample)

% 进行噪声方差已知情况下的低精度量化连续波GNOMP目标检测与估计

% 输入：
%   Tau4：等效阈值，三个维度的大小分别为：N(M)，T，4
%         其中第三个维度顺序为（实部上门限，实部下门限，虚部上门限，虚部下门限）
%   sigma：噪声方差，是一个标量
%   Phy_matrix：压缩矩阵，大小为[M,N]，非压缩场景则为单位矩阵
%   zeta：初始的信号或者时变阈值，大小为[N(M),1]
%   tau: 量化器的门限（除去最大的正无穷与最小的负无穷）
%   threshold：设置的虚警率对应的阈值
%   R_s：NOMP步骤单循环次数，默认为1
%   R_c：NOMP循环修正步骤的循环次数，默认为1
%   OverSample：过采样倍数，默认为4

% 输出：
%   OmegaList：所有目标的频率估计值，是一个列向量
%   XList：所有目标的幅度估计值（与频率一一对应），是一个列向量

if ~exist('R_s','var'), R_s = 1;
elseif isempty(R_s), R_s = 1; end

if ~exist('R_c','var'), R_c = 1;
elseif isempty(R_c), R_c = 1; end

if ~exist('OverSample','var'), OverSample = 4;
elseif isempty(OverSample), OverSample = 4; end


N = size(Phy_matrix,2);
OmegaList = [];
XList = [];
Tau4_r = Tau4;
zeta_orig = zeta;
k = 1;
while true

    % 使用RAO检测器判断是否有目标存在并且得到频率粗估计值
    omega_hat = DetectNew(Tau4_r,sigma,zeta,tau,OverSample);
    if k>K
        break
    end

    % 由频率的粗估计值得到幅度的粗估计值
    g_hat = ref_amp(Tau4_r, omega_hat, 0, sigma, Phy_matrix);

    % 对粗估计值进行牛顿修正
    [omega_est,x_est,Tau4_r] = RefineOne(Tau4_r,sigma,omega_hat,g_hat,Phy_matrix,R_s);
 
    % 将牛顿修正后的频率值和幅度值放入目标向量中
    OmegaList = [OmegaList;omega_est];
    XList = [XList;x_est];

    % 对向量中的每一个目标参数进行循环牛顿修正
    % [OmegaList,XList] = RefineAll(OmegaList,XList,Tau4_r,Phy_matrix,sigma,R_s,R_c);

    % 固定所有的信号频率，对幅度进行一次统一的修正
    % [XList,~, Tau4_r] = ref_amp(Tau4, OmegaList , XList, sigma, Phy_matrix);
    
    k = k+1;
    % 将现在估计出的所有目标对应的信号分量放入zeta中，以便下一次估计。
    zeta = zeta_orig+exp(1j*(0:N-1)'*OmegaList')*XList;
end
end