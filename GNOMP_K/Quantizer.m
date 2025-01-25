function [y_R_I_q,Tau4,y_vec_RI] = Quantizer(y_inf, Stepsize, B)
% 对信号进行量化

% 输入：
%   y_inf：未量化的原始信号
%   Stepsize：量化器的步长
%   B：量化器的比特数

% 输出：
%   y_R_I_q：量化后的信号
%   Tau4：量化后信号每个值对应的等价门限
%   y_vec_RI：使用0，1，2……代替量化后信号值处于第几个门限中

Delta_max = 2^(B-1)*Stepsize;
Delta_min = -Delta_max;
y_R1 = floor((real(y_inf)-Delta_min)/Stepsize);
y_R1(real(y_inf)>=Delta_max) = 2^B-1;
y_R1(real(y_inf)<=Delta_min) = 0;
y_I1 = floor((imag(y_inf)-Delta_min)/Stepsize);
y_I1(imag(y_inf)>=Delta_max) = 2^B-1;
y_I1(imag(y_inf)<=Delta_min) = 0;
y_vec_RI = [y_R1(:);y_I1(:)];
y_R_q1 = Delta_min + (y_R1+0.5)* Stepsize;
y_I_q1 = Delta_min + (y_I1+0.5)* Stepsize;
y_R_I_q = y_R_q1+1j*y_I_q1;

upper_th_real_g = Delta_min+Stepsize*(y_R1+1);
upper_th_real_g(y_R1==2^B-1) = inf;
lower_th_real_g = Delta_min+Stepsize*y_R1;
lower_th_real_g(y_R1==0) = -inf;

% obtain the upper and lower thresholds of the imaginary part
upper_th_imag_g = Delta_min+Stepsize*(y_I1+1);
upper_th_imag_g(y_I1==2^B-1) = inf;
lower_th_imag_g = Delta_min+Stepsize*y_I1;
lower_th_imag_g(y_I1==0) = -inf;

Tau4 = cat(3,upper_th_real_g,lower_th_real_g,upper_th_imag_g,lower_th_imag_g);

end