function [out, dbg_opt] = BiGLSE_vm_seq(y, z0, Parameter, Init, dbg_ipt)
 
%% Initilization
K_MAX = Parameter.K_MAX; P_fa = Parameter.P_fa;
N = Init.N_len;
w = (0:N-1)'/N;  w2 = w.^2; n = (0:N-1)';

Iter_in_max = Init.Iter_in_max; Iter_ref_max = Init.Iter_ref_max;
NMSE_iter_final = nan(Iter_in_max,Iter_ref_max);
omega_hat_list = nan(K_MAX, Iter_ref_max);
omega_amp_list = nan(K_MAX, Iter_ref_max);

B = Init.Bit; 
if B~=inf
    y_min = Init.y_min;
    alpha = Init.alpha;
end
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
xhat = xhat(sort_idx); omega_amp = abs(xhat); % kappa_list = kappa_list(sort_idx);
dbg_opt.omega_hat_unsel = omega_hat; dbg_opt.omega_amp_unsel = omega_amp;
mean_p = omega_hat; kappa_p = 1/(1/12*2*pi/(8*N)); % kappa_p = 1e4 % kappa_p = kappa_list; 
% figure(5)
% clf
% stem(dbg_ipt.theta, abs(dbg_ipt.x))
% hold on
% stem(omega_hat, abs(xhat));

%% algorithm iteration
L = length(omega_hat);
for iter_ot = 1:Iter_ref_max

    tau_x = 1e4*ones(L,1); xhat_damp = xhat;
    theta_hat = omega_hat; theta_hat_old = omega_hat; kappa_theta = 1e0*ones(L,1);
    mu_x0 = 0; sigvar0 = 1e3; pivar0 = 0.5*ones(L,1); sigma_w0 = 1e3;
    rhat_old = zeros(L,1); tau_r_old = 1e4*ones(L,1);
    qhat_old = omega_hat; kappa_q_old = 1e1*ones(L,1);
    shat = zeros(N,1); shat_old = zeros(N,1); tau_s_old = 1e3*ones(N,1);
    tau_p_old = 1e3*ones(N,1);
    

    for iter = 1:Iter_in_max
        
        if iter==7
            stop = 1;
        end
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
            % if B~=1
            %     v_B_ext = tau_z.*tau_p./(tau_p-tau_z);
            %     y_tilde = v_B_ext.*(zhat./tau_z-phat./tau_p);
            %     sigma_w0 = (norm(y_tilde-zhat)^2 + sum(tau_z))/N;
            % else
            sigma_w0 = Init.sigma_w;
            % end
        end

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
threshold = -sigma_w0*log(P_fa);
det = abs(omega_amp).^2;
omega_hat(det<=threshold) = []; xhat(det<=threshold) = [];
[omega_hat, sort_idx] = sort(omega_hat,'ascend'); xhat = xhat(sort_idx); 
if B==1
    NMSE_iter_final(Iter_in_max,Iter_ref_max) = norm(z0- c*exp(1j*n*omega_hat')/sqrt(N)*xhat)/norm(z0);
else
    NMSE_iter_final(Iter_in_max,Iter_ref_max) = norm(z0-exp(1j*n*omega_hat')/sqrt(N)*xhat)/norm(z0);
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