function lambda = h_func(tau, x, sigma)
% 计算h函数

% 输入：
%   tau: 量化器的门限（除去最大的正无穷与最小的负无穷）
%   x: 一个列向量
%   sigma：噪声方差

% 输出：
%   lambda: x中每个值对应的h函数，是一个列向量



% calculate the diagonal elements of the FIM
thr_upper = [tau;inf];
thr_lower = [-inf;tau];
N_x = length(x);
sigma_ext = sigma*ones(1,N_x);

t_upper = (thr_upper-x')./(sqrt(sigma_ext/2));
t_lower =(thr_lower-x')./(sqrt(sigma_ext/2));
%     Prob0 = normcdf(t_upper)-normcdf(t_lower);
Prob = 0.5 * erfc( - (1/sqrt(2))*t_upper ) - 0.5 * erfc( - (1/sqrt(2))*t_lower );
Diffpdf = (normpdf(t_upper)-normpdf(t_lower)).^2;
Prob(1,:) = 0.5 * erfc( - (1/sqrt(2))*t_upper(1,:) );
% 0.5*(1-erfc())=0.5erf()

Prob(end,:) = 0.5 * erfc( (1/sqrt(2))*t_lower(end,:) );
ratio = Diffpdf./Prob;
ratio(isinf(ratio)) = 0;
ratio(isnan(ratio)) = 0;
%     Prob0./Proba
lambda = sum(ratio);
lambda = lambda';
end
