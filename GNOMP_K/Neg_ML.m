function obj_value = Neg_ML(z,sigma,Tau4)
% 计算负的目标函数值

% 输入：
%   z：信号估计值，是一个列向量
%   sigma：噪声方差，是一个标量
%   Tau4：等效阈值，三个维度的大小分别为：N(M)，T，4
%         其中第三个维度顺序为（实部上门限，实部下门限，虚部上门限，虚部下门限）

% 输出：
%   obj_value：负的目标函数值

if ndims(Tau4)==3
    upper_th_real = Tau4(:,:,1);
    lower_th_real = Tau4(:,:,2);
    upper_th_imag = Tau4(:,:,3);
    lower_th_imag = Tau4(:,:,4);
else
    upper_th_real = Tau4(:,1);
    lower_th_real = Tau4(:,2);
    upper_th_imag = Tau4(:,3);
    lower_th_imag = Tau4(:,4);
end

upper_th_real = upper_th_real(:);
lower_th_real = lower_th_real(:);
upper_th_imag = upper_th_imag(:);
lower_th_imag = lower_th_imag(:);

z = z(:);
sigma_std = sqrt(sigma/2);

real_var_upper = (upper_th_real-real(z))/sigma_std;
R_real_upper = find(upper_th_real==inf);
real_var_lower = (lower_th_real-real(z))/sigma_std;
R_real_lower = find(lower_th_real==-inf);

imag_var_upper = (upper_th_imag-imag(z))/sigma_std;
imag_var_lower = (lower_th_imag-imag(z))/sigma_std;
I_real_upper = find(upper_th_imag==inf);
I_real_lower = find(lower_th_imag==-inf);
%     prob_real = normcdf(real_var);
prob_real = 0.5*erfc(-real_var_upper/sqrt(2))-0.5*erfc(-real_var_lower/sqrt(2));
prob_real(R_real_upper) = 1-0.5*erfc(-real_var_lower(R_real_upper)/sqrt(2));
prob_real(R_real_lower) = 0.5*erfc(-real_var_upper(R_real_lower)/sqrt(2));

prob_imag = 0.5*erfc(-imag_var_upper/sqrt(2))-0.5*erfc(-imag_var_lower/sqrt(2));
prob_imag(I_real_upper) = 1-0.5*erfc(-imag_var_lower(I_real_upper)/sqrt(2));
prob_imag(I_real_lower) = 0.5*erfc(-imag_var_upper(I_real_lower)/sqrt(2));

log_real = log(prob_real);
log_imag = log(prob_imag);

real_eq_zero = find(prob_real==0);
real_lower_abs = min(abs(real_var_lower(real_eq_zero)),abs(real_var_upper(real_eq_zero)));
% log_real(real_eq_zero) = -real_lower_abs.^2......
%     -0.5*log(2*pi)-log(real_lower_abs);
log_real(real_eq_zero) = -real_lower_abs.^2/2......
    -0.5*log(2*pi)-log(real_lower_abs);
imag_eq_zero = find(prob_imag==0);
imag_lower_abs = min(abs(imag_var_lower(imag_eq_zero)),abs(imag_var_upper(imag_eq_zero)));
% log_imag(imag_eq_zero) = -imag_lower_abs.^2......
%     -0.5*log(2*pi)-log(imag_lower_abs);
log_imag(imag_eq_zero) = -imag_lower_abs.^2/2......
    -0.5*log(2*pi)-log(imag_lower_abs);

obj_value = -sum(log_real)-sum(log_imag);
end