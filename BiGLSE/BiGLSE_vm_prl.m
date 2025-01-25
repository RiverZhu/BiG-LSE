function [out, dbg_opt] = BiGLSE_vm_prl(y, z0, a_theta, Init, dbg_ipt)
%% find initilization point
lambda_reshold = 0.1;  
Iter_max = Init.iter_max; Iter_in_max = Init.Iter_in_max; Iter_ref_max = Init.Iter_ref_max;
adaptStep = Init.adaptStep; step_default = Init.step_default;
stepMin = Init.stepMin; stepMax = Init.stepMax; stepWindow = Init.stepWindow;
stepIncr = Init.stepIncr; stepDecr = Init.stepDecr; step= Init.step_default;
N = Init.N_len; L = Init.L_len; B = Init.Bit;
if B~=inf
    y_min = Init.y_min;
    alpha = Init.alpha;
end

w = (0:N-1)'/N; w2 = w.^2; n = (0:N-1)';
theta_g = (0:2*pi/L:2*pi/L*(L-1))'; phase_rot = 1; % phase_rot = exp(-1j*(N-1)/2*theta_g);
val_x = inf; val_d = inf; down_ctn = 0; up_ctn = 0; down_flag = 0;
val_list = nan(Iter_max,1); NMSE_iter_pre = nan(Iter_max,1);

% initial algorithm parameter
sigma_w0 = 1e3;
pivar0 = 0.5; mu_x0 = 0; sigvar0 = 1e3;
sigma_d = (N/L*pi)^2/3; mean_d = 0;
xhat = zeros(L,1); tau_x = 1e3*ones(L,1);
dhat = zeros(L,1); tau_d = 1e3*ones(L,1);
shat = zeros(N,1);
xhat_damp = zeros(L,1);
shat_old = zeros(N,1); tau_s_old = 1e3*ones(N,1);
rhat_old = zeros(L,1); tau_r_old = 1e3*ones(L,1);
qhat_old = zeros(L,1); tau_q_old = 1e3*ones(L,1);
dhat_old = zeros(L,1);
tau_p_old = 1e3*ones(N,1);

%% Outer iteration 1
for iter=1:Iter_max
    
    dhat2 = dhat.^2 + eps;
    xhat2 = max(abs(xhat).^2,1e-10);
    
    tau_p_bar = mean(tau_x)*L/N+1/N*w2*(dhat2'*tau_x+xhat2'*tau_d);
    tau_p = tau_p_bar + 1/N*w2*(tau_x'*tau_d);
    tmp1 = ifft(phase_rot.*xhat);
    tmp2 = ifft(phase_rot.*xhat.*dhat);
    phat = tmp1(1:N)*L/sqrt(N)+1j*w.*tmp2(1:N)*L/sqrt(N)-shat.*tau_p_bar;
    tau_p = step*tau_p + (1-step)*tau_p_old;
    tau_p_old = tau_p;
    
    % Output channel ~ CN(0, sigma_w)
    if B==inf
        tau_z = sigma_w0.*tau_p./(sigma_w0 + tau_p);
        zhat = tau_z.*(y./sigma_w0 + phat./tau_p);
    else
        z_A_ext_real = [real(phat);imag(phat)];
        v_A_ext_real = [tau_p;tau_p]/2;
        [z_B_post_real, v_B_post_real] = GaussianMomentsComputation(y, z_A_ext_real, v_A_ext_real, y_min, B, alpha, sigma_w0/2);
        tau_z = v_B_post_real(1:N)+v_B_post_real(N+1:end); % v_B_post = mean(v_B_post)*ones(M,1);
        zhat = z_B_post_real(1:N)+1j*z_B_post_real(N+1:end);
    end
    
    if B==inf
        val_z = sum(log(2*pi*sigma_w0) + (tau_p + abs(phat-y).^2)./sigma_w0);
    else
        t_lower = y_min + alpha*y;
        t_upper = y_min + alpha*(y+1);
        t_lower(y == 0) = -inf;
        t_upper(y == 2^B-1) = +inf;
        w_up = (t_upper - z_A_ext_real)./sqrt(sigma_w0 + v_A_ext_real);
        w_dw = (t_lower - z_A_ext_real)./sqrt(sigma_w0 + v_A_ext_real);
        Prob = (erfc(w_dw) - erfc(w_up))/2;
        val_z = -log(Prob);
        val_z = sum(val_z(1:N)+val_z(N+1:end));
    end
    val = val_z + val_x + val_d;
    
    if B==inf
        sigma_w0 = (norm(y-zhat)^2 + sum(tau_z))/N;
    else
%         if B~=1
%             v_B_ext = tau_z.*tau_p./(tau_p-tau_z);
%             y_tilde = v_B_ext.*(zhat./tau_z-phat./tau_p);
%             sigma_w0 = (norm(y_tilde-zhat)^2 + sum(tau_z))/N;
%         else
            sigma_w0 = Init.sigma_w;
%         end
    end
    
    if adaptStep
        stopIdx = iter-1;
        startIdx = max(stopIdx-stepWindow,1);
        valMax = max(val_list(startIdx:stopIdx));
    end
    pass = (~adaptStep) || (iter==1) || (val <= valMax);
    if pass
        step = min(stepIncr*max(step, stepMin),stepMax);
        if adaptStep==0
            step = step_default;
        end
        down_ctn = down_ctn + 1;
        up_ctn = 0;
    else
        % stepMax = max(stepMin,maxStepDecr*stepMax);
        step = min(max(stepMin, stepDecr*step),stepMax);
        up_ctn = up_ctn + 1;
        down_ctn = 0;
        % NMSE_iter_pre(iter:end) = NMSE_iter_pre(iter-1);
    end
    if down_ctn>=3
        down_flag = 1;
        xhat_res = xhat; dhat_res = dhat;
        lambda_res = lambda; % for debug
    end
    if up_ctn<=2 && down_flag
        xhat_res = xhat; dhat_res = dhat;
        lambda_res = lambda; iter_res = iter;
    else
        down_flag = 0;
    end
    val_list(iter) = val;
    
    tau_s = (1-tau_z./tau_p)./tau_p;
    tau_s(tau_s<=0) = 1e8;
    shat = (zhat-phat)./tau_p;
    tau_s = step*tau_s + (1-step)*tau_s_old;
    shat = step*shat + (1-step)*shat_old;
    tau_s_old = tau_s;
    shat_old = shat;
    
    xhat_damp = step*xhat + (1-step)*xhat_damp;
    tau_r = 1./(mean(tau_s)+1/N*(w2'*tau_s)*dhat2);
    rhat = xhat_damp + tau_r.*(conj(phase_rot).*fft(shat,L)/sqrt(N) - 1j*a_theta'*(w.*shat).*dhat- w2'*(abs(shat).^2)/N*xhat.*tau_d);
    rhat = step*rhat+(1-step)*rhat_old;
    rhat_old = rhat;
    tau_r = step*tau_r+(1-step)*tau_r_old;
    tau_r_old = tau_r;
    
    tau_q = 1./(2*w2'*tau_s/N*(abs(xhat2).^2+eps));
    % tau_q = min(tau_q,1e8);
    qhat = dhat + 2*tau_q.*imag(conj(xhat).*( conj(phase_rot).*fft(w.*shat,L)/sqrt(N)));
    qhat = qhat.*(abs(qhat)<=pi*N/L);
    qhat = step*qhat+(1-step)*qhat_old;
    qhat_old = qhat;
    tau_q = step*tau_q+(1-step)*tau_q_old;
    tau_q_old = tau_q;
    
    % Input channel x ~ BG(lambda, mean_x, sigma_x)
    [lambda,m_t_x,V_t_x] = denoise_input_EM(pivar0,mu_x0,sigvar0,rhat,tau_r);
    lambda = max(lambda,1e-8);
    lambda = min(lambda,1-1e-8);
    xhat = lambda.*m_t_x;
    tau_x = lambda.*(abs(m_t_x).^2+V_t_x)-abs(xhat).^2+eps;
    % EM learning and KL-divergence to x
    xvar_over_sigvar = tau_x./sigvar0;
    D = -(log(xvar_over_sigvar) + 1 - xvar_over_sigvar - abs(xhat-mu_x0).^2./sigvar0);
    val_x = sum((1-lambda).*log((1-lambda)./(1-pivar0)) + lambda.*(D + log(lambda./pivar0)));
    pivar0 = lambda;
    sigvar0 = lambda'*((abs(mu_x0-m_t_x)).^2+V_t_x)/sum(lambda);
    mu_x0 = lambda'*m_t_x/sum(lambda);
    
    % % Input channel d ~ N(mean_d, sigma_d)
    tau_d = tau_q.*sigma_d./(sigma_d+tau_q);
    dhat = tau_d.*(mean_d./sigma_d+qhat./tau_q);
    dhat = step*dhat + (1-step)*dhat_old;
    dhat_old = dhat;
    % KL-divergence to d
    val_d = sum(0.5*(log(sigma_d./tau_d)+(tau_d./sigma_d-1)+(dhat-mean_d).^2./(sigma_d)));
    
    % compute NMSE/dNMSE for output
    theta_est = wrapTo2Pi(theta_g+dhat/N);
    zhat = exp(1j*n*theta_est')/sqrt(N)*xhat;
    if B==1
        debias = @(xhatvar,xvar) xhatvar'*xvar/(xhatvar'*xhatvar+eps);
        c = debias(zhat,z0);
        NMSE_iter_pre(iter) = norm(z0-zhat*c)/norm(z0);
    else
        NMSE_iter_pre(iter) = norm(z0-zhat)/norm(z0);
    end
    
end

xhat = xhat_res;
dhat = dhat_res;
lambda = lambda_res;
% export variable
out.NMSE_iter_pre = NMSE_iter_pre;

%% demension reduction
omega_ampt = abs(xhat); lambda_t = lambda; % lambda_t = ones(L,1);
for num = 1:N/3
    [~, idx] = max(omega_ampt);
    if idx==1
        omega_ampt([L,1,2])=0; lambda_t([L,2])=0;
    elseif idx==L
        omega_ampt([L-1,L,1])=0; lambda_t([L-1,1])=0;
    else
        if omega_ampt(idx-1)>omega_ampt(idx+1)
            lambda_t(idx+1)=0;
        else
            lambda_t(idx-1)=0;
        end
        omega_ampt(idx-1:idx+1)=0;
    end
end
dict_idx = find( lambda_t>lambda_reshold);
omega_hat = wrapTo2Pi( theta_g(dict_idx) + dhat(dict_idx)/N );
xhat = xhat(dict_idx);

L = length(omega_hat);
NMSE_iter_final = nan(Iter_in_max,Iter_ref_max);
omega_hat_list = nan(L, Iter_ref_max);
omega_amp_list = nan(L, Iter_ref_max);
mean_p = omega_hat; kappa_p = 1/(1/12*2*pi/(8*N)); % kappa_p = 1/(1/12*2*pi/(8*N));
% figure(5)
% clf
% stem(dbg_ipt.theta, abs(dbg_ipt.x))
% hold on
% stem(omega_hat, abs(xhat));

%% algorithm iteration
step = step_default;
for iter_ot = 1:Iter_ref_max
    
    tau_x = 1e4*ones(L,1); xhat_damp = xhat;
    theta_hat = omega_hat; theta_hat_old = omega_hat; kappa_theta = 1e0*ones(L,1);
    mu_x0 = 0; sigvar0 = 1e3; pivar0 = 0.5*ones(L,1); sigma_w0 = 1e3;
    rhat_old = zeros(L,1); tau_r_old = 1e4*ones(L,1);
    qhat_old = omega_hat; kappa_q_old = 1e1*ones(L,1);
    shat = zeros(N,1); shat_old = zeros(N,1); tau_s_old = 1e3*ones(N,1);
    tau_p_old = 1e3*ones(N,1);
    
    
    for iter = 1:Iter_in_max
        
        bsl_mat = besseli(n*ones(1,L),ones(N,1)*kappa_theta.',1)./besseli(0,ones(N,1)*kappa_theta.',1)+eps;
        bsl_mat2 = bsl_mat.^2;
        A_hat = exp(1j*n*theta_hat.').*bsl_mat/sqrt(N);
        tau_A = (1 - bsl_mat2)/N;
        xhat2 = max(abs(xhat).^2,eps); % 1e-10
        
        tau_p_bar = abs(A_hat).^2*tau_x + tau_A*abs(xhat).^2;
        tau_p = tau_p_bar + tau_A*tau_x;
        phat = A_hat*xhat - shat.*tau_p_bar;
        tau_p = step*tau_p + (1-step)*tau_p_old;
        tau_p_old = tau_p;
        
        % Output channel distribution
        if B==inf
            tau_z = sigma_w0.*tau_p./(sigma_w0 + tau_p);
            zhat = tau_z.*(y./sigma_w0 + phat./tau_p);
        else
            z_A_ext_real = [real(phat);imag(phat)];
            v_A_ext_real = [tau_p;tau_p]/2;
            [z_B_post_real, v_B_post_real] = GaussianMomentsComputation(y, z_A_ext_real, v_A_ext_real, y_min, B, alpha, sigma_w0/2);
            tau_z = v_B_post_real(1:N)+v_B_post_real(N+1:end); % v_B_post = mean(v_B_post)*ones(M,1);
            zhat = z_B_post_real(1:N)+1j*z_B_post_real(N+1:end);
        end
        if B==inf
            sigma_w0 = (norm(y-zhat)^2 + sum(tau_z))/N;
        else
            sigma_w0 = Init.sigma_w;
        end
        % sigma_w0 = Init.sigma_w;
        
        tau_s = (1-tau_z./tau_p)./tau_p;
        tau_s(tau_s<=0) = 1e8;
        shat = (zhat-phat)./tau_p;
        tau_s = step*tau_s + (1-step)*tau_s_old;
        shat = step*shat+(1-step)*shat_old;
        tau_s_old = tau_s;
        shat_old = shat;
        
        xhat_damp = step*xhat + (1-step)*xhat_damp;       
        tau_r = 1./(abs(A_hat.').^2*tau_s);
        rhat = xhat_damp + tau_r.*(A_hat'*shat - xhat.*tau_A.'*abs(shat).^2);
        rhat = step*rhat+(1-step)*rhat_old;
        rhat_old = rhat;
        tau_r = step*tau_r+(1-step)*tau_r_old;
        tau_r_old = tau_r;
        
        kappa_q = 2*N*tau_s.'*w2*(xhat2+eps);
        kappa_q = max(kappa_q, 1e8);
        qhat = theta_hat + 2*imag(conj(xhat).*( A_hat'*(w.*shat)))./kappa_q/N;
        qhat = step*qhat + (1-step)*qhat_old;
        qhat_old = qhat;
        kappa_q = step*kappa_q + (1-step)*kappa_q_old;
        kappa_q_old = kappa_q;
        
        % Input channel x ~ CN(mean_x, sigma_x)
        [lambda,m_t,V_t] = denoise_input_EM(pivar0,mu_x0,sigvar0,rhat,tau_r);
        % lambda = min(lambda,1-1e-14);
        xhat = lambda.*m_t;
        tau_x = lambda.*(abs(m_t).^2+V_t)-abs(xhat).^2;
        % EM learning to x
        pivar0 = lambda;
        sigvar0 = lambda'*((abs(mu_x0-m_t)).^2+V_t)/sum(lambda);
        mu_x0 = lambda'*m_t/sum(lambda);
        
        % Input channel d ~ N(mean_d, sigma_d)
        f_cmp = kappa_p.*exp(1j*mean_p) + kappa_q.*exp(1j*qhat);
        kappa_theta = abs(f_cmp);
        theta_hat = wrapTo2Pi(angle(f_cmp));
        theta_hat = step*theta_hat + (1-step)*theta_hat_old;
        theta_hat_old = theta_hat;
        
        % parameter learning
        mean_p = theta_hat;
        
        % theta_est = omega_hat + dhat/N;
        zhat = exp(1j*n*theta_hat')/sqrt(N)*xhat;
        if B==1
            debias = @(xhatvar,xvar) xhatvar'*xvar/(xhatvar'*xhatvar+eps);
            c = debias(zhat,z0);
            NMSE_iter_final(iter,iter_ot) = norm(z0-zhat*c)/norm(z0);
        else
            NMSE_iter_final(iter,iter_ot) = norm(z0-zhat)/norm(z0);
        end
        
    end
    
    omega_hat = theta_hat;
    omega_amp = abs(xhat);
    zhat = exp(1j*n*omega_hat')/sqrt(N)*xhat;
    NMSE_iter_final(iter,iter_ot) = norm(z0-zhat)/norm(z0);
    omega_hat_list(1:L,iter_ot) = omega_hat;
    omega_amp_list(1:L,iter_ot) = abs(xhat);
    
end

% CFAR detector
P_fa = 1e-6;
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
dbg_opt.omega_hat_aft = omega_hat_list;
dbg_opt.omega_amp_aft = omega_amp_list;
dbg_opt.xhat = xhat;
dbg_opt.threshold = threshold; dbg_opt.sigma_w0 = sigma_w0;

end