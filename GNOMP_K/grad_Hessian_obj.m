function [grad, Hessian] = grad_Hessian_obj(x_var,omega, Tau4_r, sigma, Phy_matrix)
% 计算幅度向量的梯度和hessian

% 输入：
%   x_var：幅度估计值的实部和虚部组成的列向量（[real(x);imag(x)]）
%   omega：频率估计值
%   Tau4_r：等效阈值，三个维度的大小分别为：N(M)，T，4
%           其中第三个维度顺序为（实部上门限，实部下门限，虚部上门限，虚部下门限）
%   sigma：噪声方差，是一个标量
%   Phy_matrix：压缩矩阵，大小为[M,N]，非压缩场景则为单位矩阵

% 输出：
%   grad：x_var的梯度，是一个列向量
%   Hessian：x_var的Hessian矩阵

    [M,N] = size(Phy_matrix);
    idx_all = (0:N-1)';
    A_mat = Phy_matrix*exp(1j*idx_all*omega');

    A_ext_real = [real(A_mat),-imag(A_mat)];
    A_ext_imag = [imag(A_mat),real(A_mat)];

    z_real = A_ext_real*x_var;
    z_imag = A_ext_imag*x_var;
%     z_est = z_real+1j*z_imag;
    sigma_std = sqrt(sigma/2);

    upper_th_real = Tau4_r(:,1);
    lower_th_real = Tau4_r(:,2);
    upper_th_imag = Tau4_r(:,3);
    lower_th_imag = Tau4_r(:,4);
    
    
    real_var_upper = (upper_th_real-z_real)/sigma_std;
    real_var_lower = (lower_th_real-z_real)/sigma_std;
    R_real_upper = find(upper_th_real==inf);
    R_real_lower = find(lower_th_real==-inf);

    imag_var_upper = (upper_th_imag-z_imag)/sigma_std;
    imag_var_lower = (lower_th_imag-z_imag)/sigma_std;
    I_real_upper = find(upper_th_imag==inf);
    I_real_lower = find(lower_th_imag==-inf);

    prob_real = 0.5*erfc(-real_var_upper/sqrt(2))-0.5*erfc(-real_var_lower/sqrt(2));
    prob_real(R_real_upper) = 1-0.5*erfc(-real_var_lower(R_real_upper)/sqrt(2));
    prob_real(R_real_lower) = 0.5*erfc(-real_var_upper(R_real_lower)/sqrt(2));

    prob_imag = 0.5*erfc(-imag_var_upper/sqrt(2))-0.5*erfc(-imag_var_lower/sqrt(2));
    prob_imag(I_real_upper) = 1-0.5*erfc(-imag_var_lower(I_real_upper)/sqrt(2));
    prob_imag(I_real_lower) = 0.5*erfc(-imag_var_upper(I_real_lower)/sqrt(2));
    
    real_eq_zero = find(prob_real==0);
    imag_eq_zero = find(prob_imag==0);

%     prob_real(real_eq_zero) = normpdf(real_var_lower(real_eq_zero))./real_var_lower(real_eq_zero)...
%         -normpdf(real_var_upper(real_eq_zero))./real_var_upper(real_eq_zero);
%     prob_imag(imag_eq_zero) = normpdf(imag_var_lower(imag_eq_zero))./imag_var_lower(imag_eq_zero)...
%         -normpdf(imag_var_upper(imag_eq_zero))./imag_var_upper(imag_eq_zero);
%     real_eq_zero = find(prob_real==0);
%     imag_eq_zero = find(prob_imag==0);    

    phy_real = normpdf(real_var_upper)-normpdf(real_var_lower);
    phy_image = normpdf(imag_var_upper)-normpdf(imag_var_lower);
    R_value = -1*(phy_real./prob_real);
    I_value = -1*(phy_image./prob_imag);
    R_value(real_eq_zero) = sign(real_var_upper(real_eq_zero)).*min(abs(real_var_upper(real_eq_zero)),abs(real_var_lower(real_eq_zero)));
    I_value(imag_eq_zero) = sign(imag_var_upper(imag_eq_zero)).*min(abs(imag_var_upper(imag_eq_zero)),abs(imag_var_lower(imag_eq_zero)));
    
    x_Gra_R_value = R_value.*(A_ext_real/sigma_std);
    x_Gra_I_value = I_value.*(A_ext_imag/sigma_std);
    grad = -(sum(x_Gra_R_value)+sum(x_Gra_I_value))';

    %% Hessian
    R_1_value = -R_value.^2;
    I_1_value = -I_value.^2;

    phy_d_real = normpdf(real_var_upper).*(-1*real_var_upper)-...
    normpdf(real_var_lower).*(-1*real_var_lower);
    phy_d_real(R_real_lower) = normpdf(real_var_upper(R_real_lower)).*(-1*real_var_upper(R_real_lower));
    phy_d_real(R_real_upper) = -normpdf(real_var_lower(R_real_upper)).*(-1*real_var_lower(R_real_upper));
    
    phy_d_image = normpdf(imag_var_upper).*(-1*imag_var_upper)-...
    normpdf(imag_var_lower).*(-1*imag_var_lower);
    phy_d_image(I_real_upper) = -normpdf(imag_var_lower(I_real_upper)).*(-1*imag_var_lower(I_real_upper));
    phy_d_image(I_real_lower) = normpdf(imag_var_upper(I_real_lower)).*(-1*imag_var_upper(I_real_lower));
    
    R_2_value = phy_d_real./prob_real; 
    I_2_value = phy_d_image./prob_imag;
    R_2_value(real_eq_zero) = min(abs(real_var_upper(real_eq_zero)).^2,abs(real_var_lower(real_eq_zero)).^2);
    I_2_value(imag_eq_zero) = min(abs(imag_var_upper(imag_eq_zero)).^2,abs(imag_var_lower(imag_eq_zero)).^2);
    x_Hes_R_value = A_ext_real'*((R_1_value+R_2_value).*A_ext_real);
    x_Hes_I_value = A_ext_imag'*((I_1_value+I_2_value).*A_ext_imag);
    
    Hessian = -(x_Hes_R_value+x_Hes_I_value)/(sigma/2);
    
%     if sum(isnan(grad))
%         a = 1;
%     end


end