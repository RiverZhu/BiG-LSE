function [omega_grad,omega_hessian] = Omega_grad_hessian(Tau4,omega,x_sym,sigma, Phy_matrix)
% 计算频率的一阶导数和二阶导数

% 输入：
%   Tau4：等效阈值，三个维度的大小分别为：N(M)，T，4
%           其中第三个维度顺序为（实部上门限，实部下门限，虚部上门限，虚部下门限）
%   omega：频率估计值
%   x_sym：幅度估计值乘以对称化参数
%   sigma：噪声方差，是一个标量
%   Phy_matrix：压缩矩阵，大小为[M,N]，非压缩场景则为单位矩阵

% 输出：
%   omega_grad：omega的梯度，是一个标量
%   Hessian：omega的Hessian，是一个标量

    T = size(Tau4,2);
    [M,N] = size(Phy_matrix);
%     A = exp(1j*(0:N-1)'*omega);
    upper_th_real = Tau4(:,:,1);
    lower_th_real = Tau4(:,:,2);
    upper_th_imag = Tau4(:,:,3);
    lower_th_imag = Tau4(:,:,4);
    
    upper_th_real = upper_th_real(:);
    lower_th_real = lower_th_real(:);
    upper_th_imag = upper_th_imag(:);
    lower_th_imag = lower_th_imag(:);
    
    
    n = (0:N-1)'-(N-1)/2;
    A_sym = exp(1j*n*omega);
    z = Phy_matrix*A_sym*x_sym;
    hg_sym = real(A_sym);
    ug_sym = imag(A_sym);
    
    z = z(:);
    
    omega_grad_1 = Phy_matrix*(n.*(-ug_sym*real(x_sym)-hg_sym*imag(x_sym))+1j*(n.*(-ug_sym*imag(x_sym)+hg_sym*real(x_sym))));
    omega_grad_R_1 = real(omega_grad_1);
    omega_grad_I_1 = imag(omega_grad_1);
    omega_grad_2 = Phy_matrix*(n.^2.*(-hg_sym*real(x_sym)+ug_sym*imag(x_sym))+1j*(n.^2.*(-hg_sym*imag(x_sym)-ug_sym*real(x_sym))));
    omega_grad_R_2 = real(omega_grad_2);
    omega_grad_I_2 = imag(omega_grad_2);

    sigma_std = sqrt(sigma/2);
    real_var_upper = (upper_th_real-real(z))/sigma_std;
    real_var_lower = (lower_th_real-real(z))/sigma_std;
    R_real_upper = find(upper_th_real==inf);
    R_real_lower = find(lower_th_real==-inf);

    imag_var_upper = (upper_th_imag-imag(z))/sigma_std;
    imag_var_lower = (lower_th_imag-imag(z))/sigma_std;
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

    phy_real = normpdf(real_var_upper)-normpdf(real_var_lower);
    phy_real(R_real_upper) = -normpdf(real_var_lower(R_real_upper));
    phy_real(R_real_lower) = normpdf(real_var_upper(R_real_lower) );
    phy_image = normpdf(imag_var_upper)-normpdf(imag_var_lower);
    phy_image(I_real_upper) = -normpdf(imag_var_lower(I_real_upper));
    phy_image(I_real_lower) = normpdf(imag_var_upper(I_real_lower));
    
    
    R_value = -1*(phy_real./prob_real);
    I_value = -1*(phy_image./prob_imag);
    R_value(real_eq_zero) = sign(real_var_upper(real_eq_zero)).*min(abs(real_var_upper(real_eq_zero)),abs(real_var_lower(real_eq_zero)));
    I_value(imag_eq_zero) = sign(imag_var_upper(imag_eq_zero)).*min(abs(imag_var_upper(imag_eq_zero)),abs(imag_var_lower(imag_eq_zero)));
    
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
    
    R_1_value = reshape(R_1_value,M,T);
    R_2_value = reshape(R_2_value,M,T);
    I_1_value = reshape(I_1_value,M,T);
    I_2_value = reshape(I_2_value,M,T);
    
    R_value = reshape(R_value,M,T);
    I_value = reshape(I_value,M,T);
    
    omega_Hes_R_value = (R_1_value+R_2_value).*omega_grad_R_1.^2+R_value.*omega_grad_R_2;
    omega_Hes_I_value = (I_1_value+I_2_value).*omega_grad_I_1.^2+I_value.*omega_grad_I_2;

    omega_hessian = -(sum(omega_Hes_R_value,'all')+sum(omega_Hes_I_value,'all'))/(sigma/2);
    
    omega_Gra_R_value = R_value.*(omega_grad_R_1/sigma_std);
    omega_Gra_I_value = I_value.*(omega_grad_I_1/sigma_std);

    omega_grad = -(sum(omega_Gra_R_value,'all')+sum(omega_Gra_I_value,'all'));


end
