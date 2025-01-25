function [omega_hat,Test_RAO_DFT] = DetectNew(Tau4_r,sigma,zeta,tau,OverSample)
% 计算噪声方差已知情况下的RAO TEST

% 输入：
%   Tau4_r：等效阈值，三个维度的大小分别为：N(M)，T，4
%           其中第三个维度顺序为（实部上门限，实部下门限，虚部上门限，虚部下门限）
%   sigma：噪声方差，是一个标量
%   zeta：初始的信号或者时变阈值，大小为[N(M),1]
%   tau: 量化器的门限（除去最大的正无穷与最小的负无穷）
%   OverSample：过采样倍数

% 输出：
%   omega_hat：过采样的RAO检测统计量中最大值对应的频率，是一个标量
%   Test_RAO_DFT：DFT格点上的RAO检测统计量的值，是一个列向量

N = size(Tau4_r,1);
N_OS = N*OverSample;
Omega_OS_Grid = (0:N_OS-1)*2*pi/N_OS;
H_1 = (h_func(tau,real(zeta),sigma)+h_func(tau,imag(zeta),sigma))/2;
H_2 = (h_func(tau,real(zeta),sigma)-h_func(tau,imag(zeta),sigma))/2;

real_upper_x = Tau4_r(:,1,1)/sqrt(sigma/2);
real_lower_x = Tau4_r(:,1,2)/sqrt(sigma/2);
imag_upper_x = Tau4_r(:,1,3)/sqrt(sigma/2);
imag_lower_x = Tau4_r(:,1,4)/sqrt(sigma/2);

real_upper_inf = real_upper_x==inf;
real_lower_inf = real_lower_x==-inf;
imag_upper_inf = imag_upper_x==inf;
imag_lower_inf = imag_lower_x==-inf;

upper_real = normpdf(real_upper_x)-normpdf(real_lower_x);
upper_imag = normpdf(imag_upper_x)-normpdf(imag_lower_x);
lower_real = 0.5*erfc(-real_upper_x/sqrt(2))-0.5*erfc(-real_lower_x/sqrt(2));
lower_real(real_upper_inf) = 1-0.5*erfc(-real_lower_x(real_upper_inf)/sqrt(2));
lower_real(real_lower_inf) = 0.5*erfc(-real_upper_x(real_lower_inf)/sqrt(2));
lower_imag = 0.5*erfc(-imag_upper_x/sqrt(2))-0.5*erfc(-imag_lower_x/sqrt(2));
lower_imag(imag_upper_inf) = 1-0.5*erfc(-imag_lower_x(imag_upper_inf)/sqrt(2));
lower_imag(imag_lower_inf) = 0.5*erfc(-imag_upper_x(imag_lower_inf)/sqrt(2));
lower_real_zero = find(lower_real==0);
lower_imag_zero = find(lower_imag==0);
phi_real = -upper_real./lower_real;
phi_imag = -upper_imag./lower_imag;
phi_real(lower_real_zero)=sign(real_upper_x(lower_real_zero)).*min(abs(real_upper_x(lower_real_zero)),abs(real_lower_x(lower_real_zero)));
phi_imag(lower_imag_zero)=sign(imag_upper_x(lower_imag_zero)).*min(abs(imag_upper_x(lower_imag_zero)),abs(imag_lower_x(lower_imag_zero)));
phi = phi_real+1j*phi_imag;

if size(H_1,1) == 1
    H_1_sum = H_1*N;
else
    H_1_sum = sum(H_1);
end
A_H_phy = fft(phi,N_OS);

H_2_ext = zeros(2*N,1);
H_2_ext(1:2:end) = H_2;
AT_H_AT = ifft(H_2_ext,N_OS)*N_OS;
Test_Rao_all = (H_1_sum.*abs(A_H_phy).^2-real(AT_H_AT.*A_H_phy.^2))...
    ./(H_1_sum.^2-abs(AT_H_AT).^2);
Test_RAO_DFT = Test_Rao_all(1:OverSample:end);
[~,idx_max] = max(Test_Rao_all);
omega_hat = Omega_OS_Grid(idx_max);
end

