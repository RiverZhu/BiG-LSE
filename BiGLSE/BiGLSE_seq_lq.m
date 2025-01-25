function [out, dbg_opt] = BiGLSE_seq_lq(y, z0, Parameter, Init, dbg_ipt)

%% Initilization
OverSampling = Parameter.OverSampling; K_MAX = Parameter.K_MAX; P_fa = Parameter.P_fa;
N = Init.N_len; L = N*OverSampling;
fac = 2; van_th = 1e-2; omega_th = 0.4*2*pi/N;
sigma_d = (fac*N/L*pi)^2/3; mean_d = 0;
w = ((0:N-1)-(N-1)/2)'/N;  w2 = w.^2;
n = (-(N-1)/2:(N-1)/2)';

Iter_in_max = Init.Iter_in_max; Iter_ref_max = Init.Iter_ref_max;
NMSE_iter_final = nan(Iter_in_max,Iter_ref_max);
omgea_hat_list = nan(K_MAX, Iter_ref_max);
omgea_amp_list = nan(K_MAX, Iter_ref_max);

B = Init.Bit;
step = Init.step_default;

%% find initilization point
if B==inf
    [omega_hat, xhat, ~] = NOMPinit(y, eye(N), K_MAX);
else
    sigma_w = Init.sigma_w; % used for initialization
    y_min = Init.y_min;
    alpha = Init.alpha;
    upper_th_g = y_min + alpha*(y+1);
    upper_th_g(y==2^B-1) = inf;
    lower_th_g = y_min + alpha*y;
    lower_th_g(y==0) = -inf;
    Tau4 = cat(3,upper_th_g(1:N),lower_th_g(1:N),upper_th_g(N+1:2*N),lower_th_g(N+1:2*N));
    zeta = 0; % zeta = Init.zeta;
    Tau4_r = Tau4-cat(3,real(zeta),real(zeta),imag(zeta),imag(zeta));
    tau = -2^(B-1)*alpha+(0:2^B)'*alpha; tau = tau(2:end-1);
    [omega_hat,xhat] = GNOMP_K(Tau4_r, sigma_w, eye(N), zeta, tau, K_MAX);
end

[omega_hat, sort_idx] = sort(omega_hat,'ascend');
xhat = xhat(sort_idx); omega_amp = abs(xhat);
dbg_opt.omega_hat_unsel = omega_hat; dbg_opt.omega_amp_unsel = omega_amp;
omega_dif = diff(omega_hat); omega_dif_idx = find( omega_dif<=(2*pi/OverSampling/N+eps) );
del_idx = omega_dif_idx(omega_amp(omega_dif_idx)<omega_amp(omega_dif_idx+1));
omega_hat(del_idx) = []; xhat(del_idx) = []; omega_amp(del_idx) = [];
dbg_opt.omega_hat_pre = omega_hat; dbg_opt.omega_amp_pre = omega_amp;

% figure(5)
% clf
% stem(dbg_ipt.theta, abs(dbg_ipt.x))
% hold on
% stem(omega_hat, abs(xhat));

%% 算法迭代
L = length(omega_hat);
for iter_ot = 1:Iter_ref_max
    
    if iter_ot~=1
        dict_idx = find(abs(xhat)<=van_th);
        omega_hat(dict_idx) = []; xhat(dict_idx) = []; omega_amp(dict_idx) = [];
        [omega_amp_list, order] = sort(omega_amp); xhat_list = xhat(order); omega_hat_list = omega_hat(order);
        arr = omega_hat_list(end); arr_xhat = xhat_list(end);
        omega_hat_list(end) = []; xhat_list(end) = []; omega_amp_list(end) = [];
        for i=1:length(omega_hat)-1
            omega_sel = omega_hat_list(end);
            if i==11
                11;
            end
            if ~sum(abs(wrapToPi(arr-omega_sel)) <= omega_th)
                arr = [arr; omega_sel];
                arr_xhat = [arr_xhat; xhat_list(end)];
            else
                idx_fus = find(abs(wrapToPi(arr-omega_sel))<=omega_th);
                if length(idx_fus)~=1
                    [~,idx_one] = min(arr(idx_fus)-omega_sel);
                    idx_fus = idx_fus(idx_one);
                end
                if ~(abs(arr_xhat(idx_fus))>=0.4 && abs(xhat_list(end))>=0.4 && abs(wrapToPi(arr(idx_fus)-omega_sel))>=0.3*2*pi/N)
                    % fre_wei = abs(arr_xhat(idx_fus)) / (abs(xhat_list(end)) + abs(arr_xhat(idx_fus)));
                    fre_wei = 1;
                    arr(idx_fus) = fre_wei*arr(idx_fus) + (1-fre_wei)*omega_sel;
                    arr_xhat(idx_fus) = arr_xhat(idx_fus) + xhat_list(end);
                else
                    arr = [arr; omega_sel];
                    arr_xhat = [arr_xhat; xhat_list(end)];
                end
            end
            omega_hat_list(end) = []; xhat_list(end) = []; omega_amp_list(end) = [];
        end
        omega_hat = arr; L = length(arr); xhat = zeros(L,1);
    end
    
    tau_x = 1e4*ones(L,1); xhat_damp = xhat;
    dhat = zeros(L,1); tau_d = 1e3*ones(L,1);
    a_theta = exp(1j*n*omega_hat')/sqrt(N);
    mu_x0 = 0; sigvar0 = 1e4; pivar0 = 0.5*ones(L,1); sigma_w0 = 1e3;
    rhat_old = zeros(L,1); tau_r_old = 1e2*ones(L,1);
    qhat_old = zeros(L,1); tau_q_old = 1e2*ones(L,1);
    shat = zeros(N,1); shat_old = zeros(N,1); tau_s_old = 1e4*ones(N,1);
    tau_p_old = 1e4*ones(N,1);
    dhat_old = zeros(L,1);
    
    for iter = 1:Iter_in_max
        
        dhat2 = (dhat).^2  + eps;
        xhat2 = max(abs(xhat).^2,eps); % 1e-10
        
        tau_p_bar = mean(tau_x)*L/N + 1/N*w2*(dhat2'*tau_x + xhat2'*tau_d);
        tau_p =  tau_p_bar + 1/N*w2*(tau_x'*tau_d);
        phat = a_theta*xhat + 1j*w.*(a_theta*(xhat.*dhat)) - shat.*tau_p_bar;
        tau_p = step*tau_p + (1-step)*tau_p_old;
        tau_p_old = tau_p;
        
        % Output channel distribution
        if B==inf
            tau_z = sigma_w0.*tau_p./(sigma_w0 + tau_p);
            zhat = tau_z.*(y./sigma_w0 + phat./tau_p);
        else
            z_A_ext_real = [real(phat);imag(phat)];
            v_A_ext_real = [tau_p;tau_p]/2;
            zeta = 0; % zeta = [real(Init.zeta);imag(Init.zeta)];
            [z_B_post_real, v_B_post_real] = GaussianMomentsComputation_zeta(y, z_A_ext_real, v_A_ext_real, y_min, B, alpha, sigma_w0/2, zeta);
            tau_z = v_B_post_real(1:N)+v_B_post_real(N+1:end); % v_B_post = mean(v_B_post)*ones(M,1);
            zhat = z_B_post_real(1:N)+1j*z_B_post_real(N+1:end);
        end
        
        if B==inf
            sigma_w0 = (norm(y-zhat)^2 + sum(tau_z))/N;
        else
            sigma_w0 = Init.sigma_w;
        end
        
        tau_s = (1-tau_z./tau_p)./tau_p;
        tau_s(tau_s<=0) = 1e8;
        shat = (zhat-phat)./tau_p;
        tau_s = step*tau_s + (1-step)*tau_s_old;
        shat = step*shat+(1-step)*shat_old;
        tau_s_old = tau_s;
        shat_old = shat;
        
        xhat_damp = step*xhat + (1-step)*xhat_damp;
        tau_r = 1./(mean(tau_s) + (w2'*tau_s)/N*dhat2);
        rhat = xhat_damp + tau_r.*(a_theta'*shat -1j*a_theta'*(w.*shat).*dhat- w2'*(abs(shat).^2)/N*xhat.*tau_d);
        rhat = step*rhat+(1-step)*rhat_old;
        rhat_old = rhat;
        tau_r = step*tau_r+(1-step)*tau_r_old;
        tau_r_old = tau_r;
        
        tau_q = 1./(2*w2'*tau_s/N*(xhat2+eps));
        qhat = dhat + 2*tau_q.*imag(conj(xhat).*( a_theta'*(w.*shat)));
        qhat = qhat.*(abs(qhat)<=fac*pi/OverSampling);
        % qhat = qhat.*(abs(qhat)<=pi*N/L);
        
        qhat = step*qhat + (1-step)*qhat_old;
        qhat_old = qhat;
        tau_q = step*tau_q + (1-step)*tau_q_old;
        tau_q_old = tau_q;
        
        % Input channel x ~ CN(mean_x, sigma_x)
        [lambda,m_t,V_t] = denoise_input_EM(pivar0,mu_x0,sigvar0,rhat,tau_r);
        % lambda = max(lambda,1e-14);
        lambda = min(lambda,1-1e-14);
        xhat = lambda.*m_t;
        tau_x = lambda.*(abs(m_t).^2+V_t)-abs(xhat).^2;
        % EM learning to x
        pivar0 = lambda;
        sigvar0 = lambda'*((abs(mu_x0-m_t)).^2+V_t)/sum(lambda);
        mu_x0 = lambda'*m_t/sum(lambda);
        
        % Input channel d ~ N(mean_d, sigma_d)
        tau_d = tau_q.*sigma_d./(sigma_d+tau_q);
        dhat = tau_d.*(mean_d./sigma_d+qhat./tau_q);
        % [dhat_dy,tau_d_dy] = calculate_n_k(qhat_dy,pi*N/L,tau_q_dy);
        dhat = step*dhat + (1-step)*dhat_old;
        dhat_old = dhat;
        
        theta_est = omega_hat + dhat/N;
        zhat = exp(1j*n*theta_est')/sqrt(N)*xhat;
        if B==1
            debias = @(xhatvar,xvar) xhatvar'*xvar/(xhatvar'*xhatvar+eps);
            c = debias(zhat,z0);
            NMSE_iter_final(iter,iter_ot) = norm(z0-zhat*c)/norm(z0);
        else
            NMSE_iter_final(iter,iter_ot) = norm(z0-zhat)/norm(z0);
        end
        
    end
    
    omega_hat = omega_hat + dhat/N;
    omega_amp = abs(xhat);
    zhat = exp(1j*n*omega_hat')/sqrt(N)*xhat;
    if B==1
        c = debias(zhat,z0);
        NMSE_iter_final(iter,iter_ot) = norm(z0-zhat*c)/norm(z0);
    else
        NMSE_iter_final(iter,iter_ot) = norm(z0-zhat)/norm(z0);
    end
    omgea_hat_list(1:L,iter_ot) = omega_hat;
    omgea_amp_list(1:L,iter_ot) = abs(xhat);

end

% CFAR detector
threshold = -sigma_w0*log(P_fa);
det = abs(omega_amp).^2;
omega_hat(det<=threshold) = []; xhat(det<=threshold) = [];
[omega_hat, sort_idx] = sort(omega_hat,'ascend'); xhat = xhat(sort_idx);
zhat = exp(1j*n*omega_hat')/sqrt(N)*xhat;
if B==1
    NMSE_iter_final(Iter_in_max,Iter_ref_max) = norm(z0-c*zhat)/norm(z0);
else
    NMSE_iter_final(Iter_in_max,Iter_ref_max) = norm(z0-zhat)/norm(z0);
end

% output and debug
out.omega_hat = omega_hat;
out.omega_amp = xhat;
if B==1
    out.omega_amp = xhat*c;
end
out.NMSE_iter_final = NMSE_iter_final;
dbg_opt.omgea_hat_aft = omgea_hat_list;
dbg_opt.omgea_amp_aft = omgea_amp_list;
dbg_opt.xhat = xhat;
dbg_opt.threshold = threshold; dbg_opt.sigma_w0 = sigma_w0;

end