function omega_hat_refine = ref_freq(Tau4_r, omega_hat, ghat, sigma,Phy_matrix)
% 对目标的频率进行牛顿修正

% 输入：
%   Tau4_r_T：等效阈值，三个维度的大小分别为：N(M)，T，4
%           其中第三个维度顺序为（实部上门限，实部下门限，虚部上门限，虚部下门限）
%   omega_hat：频率估计值，是一个标量
%   ghat：幅度估计值，是一个标量
%   sigma：噪声方差，是一个标量
%   Phy_matrix：压缩矩阵，大小为[M,N]，非压缩场景则为单位矩阵

% 输出：
%   omega_hat_refine：频率估计值修正后的值


[M,N] = size(Phy_matrix);
T = size(Tau4_r,2);
omega0 = omega_hat;
ghat_sym = ghat*exp(1j*(N-1)/2*omega_hat);
step_size = 1;
for iter2 = 1:2

    [omega_grad,omega_hess] = Omega_grad_hessian(Tau4_r,omega0,ghat_sym,sigma, Phy_matrix);
    if omega_hess==0
        omega0 = omega0-step_size*sign(omega_grad)*pi/N/32;
    else
        Delta_omega = -omega_hess\omega_grad;
        if abs(Delta_omega)>pi/N/4
            Delta_omega = -sign(omega_grad)*pi/N/32;
            %             num = num+1
        end
        omega0 = omega0+step_size*Delta_omega;
    end
end
omega_hat_refine = omega0;
% omega_hat_refine = wrapTo2Pi(omega0) ;
end