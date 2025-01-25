function [ghat,neg_log_refine, Tau4_r_detect] = ref_amp(Tau4_r_T, omega_hat, ghat_init, sigma, Phy_matrix)
% 对目标的幅度进行牛顿修正

% 输入：
%   Tau4_r_T：等效阈值，三个维度的大小分别为：N(M)，T，4
%           其中第三个维度顺序为（实部上门限，实部下门限，虚部上门限，虚部下门限）
%   omega_hat：频率估计值，是一个标量
%   ghat_init：幅度估计值，是一个标量
%   sigma：噪声方差，是一个标量
%   Phy_matrix：压缩矩阵，大小为[M,N]，非压缩场景则为单位矩阵

% 输出：
%   ghat：幅度估计值修正后的值
%   neg_log_refine：此时对应的负的目标函数值的大小，是一个标量
%   Tau4_r_detect：修正后的等效阈值

[M,N] = size(Phy_matrix);
[K,T] = size(ghat_init);
ghat = zeros(K,T);

x_cplx = @(x_real) x_real(1:end/2)+1j*x_real(end/2+1:end);
A = Phy_matrix*exp(1j*(0:N-1)'*omega_hat');
alpha = 0.2;
beta = 0.5;
tol = 1e-6;
delta_ML = 0.01;
step_size_grad_dec = 1e-4;

for iter = 1:T
    Tau4_r = Tau4_r_T(:,iter,:);

    % 将幅度的实部和虚部放在同一个列向量中
    x0 = [real(ghat_init(:,iter));imag(ghat_init(:,iter))];
    num = 0;
    while true
        num = num+1;
        
        % 计算x0对应的梯度和Hessian，并进行牛顿修正。
        [grad, Hessian] = grad_Hessian_obj(x0,omega_hat, Tau4_r, sigma, Phy_matrix);
        z0 = A*x_cplx(x0);
        obj_current = Neg_ML(z0,sigma,Tau4_r); 
        step_size = 1;
        % 进行一步异常值保护措施
        if rank(Hessian)==size(Hessian,1)
            Deltag = -Hessian\grad;
            lambda = grad'/Hessian*grad;
            x1 = x0+step_size*Deltag;
            z1 = A*x_cplx(x1);
            obj_new = Neg_ML(z1,sigma,Tau4_r); 
            if lambda<tol
                break;
            end
            while obj_new-obj_current>alpha*step_size*grad'*Deltag
                step_size = beta*step_size;
                x1 = x0+step_size*Deltag;
                z1 = A*x_cplx(x1);
                obj_new = Neg_ML(z1,sigma,Tau4_r);
            end
            x0 = x1;
        else
            x1 = x0-step_size_grad_dec*grad;
            z1 = A*x_cplx(x1);
            obj_new = Neg_ML(z1,sigma,Tau4_r);
            while obj_new-obj_current>0
                step_size_grad_dec = beta*step_size_grad_dec;
                x1 = x0-step_size_grad_dec*grad;
                z1 = A*x_cplx(x1);
                obj_new = Neg_ML(z1,sigma,Tau4_r);
            end
            x0 = x1; 
        end
        if abs(obj_new-obj_current)<delta_ML
            break;
        end

        
    end
    ghat(:,iter) = x_cplx(x0);
    z_refine = A*ghat;
    % 生成新的目标函数值和等价阈值
    neg_log_refine = Neg_ML(z_refine,sigma,Tau4_r_T);
    Tau4_r_detect = Tau4_r_T -cat(3,real(z_refine),real(z_refine),imag(z_refine),imag(z_refine)) ;
end