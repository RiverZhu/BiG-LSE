function [out, dbg_opt] = BiGLSE_prl(y, z0, a_theta, Init, dbg_ipt)

    %% initilization parameter
    % initilize default parameter
    Iter_max = Init.iter_max; Iter_in_max = Init.Iter_in_max; Iter_ref_max = Init.Iter_ref_max; 
    adaptStep = Init.adaptStep; step_default = Init.step_default; 
    stepMin = Init.stepMin; stepMax = Init.stepMax; stepWindow = Init.stepWindow;
    stepIncr = Init.stepIncr; stepDecr = Init.stepDecr; step= Init.step_default;
    N = Init.N_len; L = Init.L_len; B = Init.Bit;
    if B~=inf
        y_min = Init.y_min;
        alpha = Init.alpha;
    end
    debug = Init.debug;
    % omega_threshold = 0.355*2*pi/N; lambda_reshold = 0.1;
    omega_threshold = 0.4*2*pi/N; lambda_reshold = 0.1;
    
    w = ((0:N-1)-(N-1)/2)'/N;  w2 = w.^2;
    theta_g = (0:2*pi/L:2*pi/L*(L-1))'; phase_rot = exp(-1j*(N-1)/2*theta_g);
    val_x = inf; val_d = inf; 
    down_ctn = 0; up_ctn = 0; down_flag = 0;
    val_list = nan(Iter_max,1);
    NMSE_iter_pre = nan(Iter_max,1);
    
    % initial algorithm parameter
    sigma_w0 = 1e2;
    pivar0 = 0.4; mu_x0 = 0; sigvar0 = 1e2;
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
        
        if iter==24
            1; % for debug
        end
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
            if B~=1
                v_B_ext = tau_z.*tau_p./(tau_p-tau_z);
                y_tilde = v_B_ext.*(zhat./tau_z-phat./tau_p);
                sigma_w0 = (norm(y_tilde-zhat)^2 + sum(tau_z))/N;
            else
                sigma_w0 = Init.sigma_w;
            end
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
        end
        if down_ctn>=3 
            down_flag = 1;
            xhat_res = xhat; dhat_res = dhat;
            lambda_res = lambda; iter_res = iter; % for debug
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
        % [dhat_t,tau_d_t] = calculate_n_k(qhat,pi*N/L,tau_q);
        dhat = step*dhat + (1-step)*dhat_old;
        dhat_old = dhat;
        % KL-divergence to d
        val_d = sum(0.5*(log(sigma_d./tau_d)+(tau_d./sigma_d-1)+(dhat-mean_d).^2./(sigma_d)));

        % compute NMSE/dNMSE for output
        theta_est = wrapTo2Pi(theta_g+dhat/N);
        n = (-(N-1)/2:(N-1)/2)';
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
    if debug==1
        dbg_opt.omega_hat_unsel = theta_est;
        dbg_opt.omega_amp_unsel = abs(xhat);
        dbg_opt.sigma_w0 = sigma_w0;
        dbg_opt.val_list = val_list;
        dbg_opt.iter_res = iter_res;
    end
    
    
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
    omega_amp = abs(xhat(dict_idx));
    L_dy = length(omega_hat); xhat_dy = zeros(L_dy,1); % xhat = randn(L_dy,1)+ 1j*randn(L_dy,1);
    NMSE_iter_final = zeros(Iter_in_max,Iter_ref_max);
    
    % For debug
    if debug
        dbg_opt.omega_hat_pre = omega_hat;
        if B==1
            dbg_opt.omega_amp_pre = abs(xhat(dict_idx)*c);
        else
            dbg_opt.omega_amp_pre = omega_amp;
        end
        dbg_opt.omgea_hat_aft = nan(L_dy, Iter_ref_max);
        dbg_opt.omgea_amp_aft = nan(L_dy, Iter_ref_max);
    end
    
    
    %% Outer iteration 2->end
    damp_fac_dy = 0.8;
    xhat_dy_damp = xhat_dy;
    
    for iter_ref = 1:Iter_ref_max
        
        if iter_ref>=3
            if iter_ref~=1
                dict_idx = find(abs(xhat_dy)<=5e-2);
                omega_hat(dict_idx) = []; xhat_dy(dict_idx) = []; omega_amp(dict_idx) = [];
            end
            
            omega_hat_list = omega_hat;
            omega_amp_list = omega_amp;
            xhat_list = xhat_dy;
            [~, idx] = max(omega_amp_list);
            arr = omega_hat_list(idx); arr_xhat = xhat_list(idx); % arr_amp = omega_amp_list(idx);
            omega_hat_list(idx) = []; xhat_list(idx) = []; omega_amp_list(idx) = [];
            for i=1:length(omega_hat)-1
                [~, idx] = max(omega_amp_list);
                omega_sel = omega_hat_list(idx);
                if ~sum(abs(wrapToPi(arr-omega_sel)) <= omega_threshold)
                    arr = [arr; omega_sel];
                    arr_xhat = [arr_xhat; xhat_list(idx)];

                else
                    idx_fus = find(abs(wrapToPi(arr-omega_sel))<=omega_threshold);
                    arr_xhat(idx_fus) = arr_xhat(idx_fus) + xhat_list(idx);
                end
                omega_hat_list(idx) = []; xhat_list(idx) = []; omega_amp_list(idx) = [];
            end
            
            omega_hat = arr; L_dy = length(arr);
            % xhat_dy = arr_xhat; xhat_dy_damp = arr_xhat;
            xhat_dy = zeros(L_dy,1); xhat_dy_damp = xhat_dy; 
        end
        
        pivar0_dy = 0.5*ones(L_dy,1);
        a_theta_dy = exp(1j*n*omega_hat')/sqrt(N);
        tau_p_dy_old = 1e2*ones(N,1);
        tau_s_dy_old = 1e2*ones(N,1);
        dhat_dy_old = zeros(L_dy,1);
        tau_x_dy = 1e2*ones(L_dy,1);
        dhat_dy = zeros(L_dy,1);
        tau_d_dy = 1e2*ones(L_dy,1);
        rhat_dy_old = zeros(L_dy,1);
        tau_r_dy_old = 1e2*ones(L_dy,1);
        qhat_dy_old = zeros(L_dy,1);
        tau_q_dy_old = 1e2*ones(L_dy,1);
        shat_dy = zeros(N,1);
        shat_dy_old = zeros(N,1);
        sigma_w0 = 1e2;

        for iter = 1:Iter_in_max

            dhat2 = (dhat_dy).^2  + eps;
            xhat2 = max(abs(xhat_dy).^2,1e-10);

            tau_p_dy_bar = mean(tau_x_dy)*L_dy/N + 1/N*w2*(dhat2'*tau_x_dy + xhat2'*tau_d_dy);
            tau_p_dy =  tau_p_dy_bar + 1/N*w2*(tau_x_dy'*tau_d_dy);
            phat_dy = a_theta_dy*xhat_dy + 1j*w.*(a_theta_dy*(xhat_dy.*dhat_dy)) - shat_dy.*tau_p_dy_bar;
            tau_p_dy = damp_fac_dy*tau_p_dy + (1-damp_fac_dy)*tau_p_dy_old;
            tau_p_dy_old = tau_p_dy;

            % Output channel distribution
            if B==inf
                tau_z_dy = sigma_w0.*tau_p_dy./(sigma_w0 + tau_p_dy);
                zhat_dy = tau_z_dy.*(y./sigma_w0 + phat_dy./tau_p_dy);
            else
                z_A_ext_real = [real(phat_dy);imag(phat_dy)];
                v_A_ext_real = [tau_p_dy;tau_p_dy]/2;
                [z_B_post_real, v_B_post_real] = GaussianMomentsComputation(y, z_A_ext_real, v_A_ext_real, y_min, B, alpha, sigma_w0/2);
                tau_z_dy = v_B_post_real(1:N)+v_B_post_real(N+1:end); % v_B_post = mean(v_B_post)*ones(M,1);
                zhat_dy = z_B_post_real(1:N)+1j*z_B_post_real(N+1:end);
            end

            if B==inf
                sigma_w0 = (norm(y-zhat_dy)^2 + sum(tau_z_dy))/N;
            else
                if B~=1
                    v_B_ext = tau_z_dy.*tau_p_dy./(tau_p_dy-tau_z_dy);
                    y_tilde = v_B_ext.*(zhat_dy./tau_z_dy-phat_dy./tau_p_dy);
                    sigma_w0 = (norm(y_tilde-zhat_dy)^2 + sum(tau_z_dy))/N;
                else
                    sigma_w0 = Init.sigma_w;
                end
            end

            tau_s_dy = (1-tau_z_dy./tau_p_dy)./tau_p_dy;
            tau_s_dy(tau_s_dy<=0) = 1e8;
            shat_dy = (zhat_dy-phat_dy)./tau_p_dy;
            tau_s_dy = damp_fac_dy*tau_s_dy + (1-damp_fac_dy)*tau_s_dy_old;
            shat_dy = damp_fac_dy*shat_dy+(1-damp_fac_dy)*shat_dy_old;
            tau_s_dy_old = tau_s_dy;
            shat_dy_old = shat_dy;

            xhat_dy_damp = damp_fac_dy*xhat_dy + (1-damp_fac_dy)*xhat_dy_damp;
            tau_r_dy = 1./(mean(tau_s_dy) + (w2'*tau_s_dy)/N*dhat2);
            rhat_dy = xhat_dy_damp + tau_r_dy.*(a_theta_dy'*shat_dy -1j*a_theta_dy'*(w.*shat_dy).*dhat_dy- w2'*(abs(shat_dy).^2)/N*xhat_dy.*tau_d_dy);
            rhat_dy = damp_fac_dy*rhat_dy+(1-damp_fac_dy)*rhat_dy_old;
            rhat_dy_old = rhat_dy;
            tau_r_dy = damp_fac_dy*tau_r_dy+(1-damp_fac_dy)*tau_r_dy_old;
            tau_r_dy_old = tau_r_dy;

            tau_q_dy = 1./(2*w2'*tau_s_dy/N*(xhat2+eps));
            qhat_dy = dhat_dy + 2*tau_q_dy.*imag(conj(xhat_dy).*( a_theta_dy'*(w.*shat_dy)));
            qhat_dy = qhat_dy.*(abs(qhat_dy)<=pi*N/L);

            qhat_dy = damp_fac_dy*qhat_dy + (1-damp_fac_dy)*qhat_dy_old;
            qhat_dy_old = qhat_dy;
            tau_q_dy = damp_fac_dy*tau_q_dy + (1-damp_fac_dy)*tau_q_dy_old;
            tau_q_dy_old = tau_q_dy;

            % Input channel x ~ CN(mean_x, sigma_x)
            [lambda,m_t,V_t] = denoise_input_EM(pivar0_dy,mu_x0,sigvar0,rhat_dy,tau_r_dy);
            lambda = max(lambda,1e-8);
            lambda = min(lambda,1-1e-8);
            xhat_dy = lambda.*m_t;
            tau_x_dy = lambda.*(abs(m_t).^2+V_t)-abs(xhat_dy).^2+eps;
            % EM learning to x
            pivar0_dy = lambda;
            sigvar0 = lambda'*((abs(mu_x0-m_t)).^2+V_t)/sum(lambda);
            mu_x0 = lambda'*m_t/sum(lambda);

            % Input channel d ~ N(mean_d, sigma_d)
            tau_d_dy = tau_q_dy.*sigma_d./(sigma_d+tau_q_dy);
            dhat_dy = tau_d_dy.*(mean_d./sigma_d+qhat_dy./tau_q_dy);
            % [dhat_dy,tau_d_dy] = calculate_n_k(qhat_dy,pi*N/L,tau_q_dy);
            dhat_dy = damp_fac_dy*dhat_dy + (1-damp_fac_dy)*dhat_dy_old;
            dhat_dy_old = dhat_dy;

            theta_est = omega_hat + dhat_dy/N;          
            zhat_dy = exp(1j*n*theta_est')/sqrt(N)*xhat_dy;
            out.zhat_dy = zhat_dy;
            if B==1
                debias = @(xhatvar,xvar) xhatvar'*xvar/(xhatvar'*xhatvar+eps);
                c = debias(zhat_dy,z0);
                NMSE_iter_final(iter,iter_ref) = norm(z0-zhat_dy*c)/norm(z0);
            else
                NMSE_iter_final(iter,iter_ref) = norm(z0-zhat_dy)/norm(z0);
            end

        end

        omega_hat = wrapTo2Pi(omega_hat + dhat_dy/N); omega_amp = abs(xhat_dy);
        % for debug
        dbg_opt.omgea_hat_aft(1:L_dy,iter_ref) = omega_hat;
        dbg_opt.omgea_amp_aft(1:L_dy,iter_ref) = omega_amp;
        
    end

    % CFAR detector
    P_fa = 1e-6;
    threshold = -sigma_w0*log(P_fa);
    det = abs(omega_amp).^2;
    omega_hat(det<=threshold) = []; xhat_dy(det<=threshold) = [];
    [omega_hat, sort_idx] = sort(omega_hat,'ascend');
    xhat_dy = xhat_dy(sort_idx); omega_amp = xhat_dy;
    zhat_dy = exp(1j*n*omega_hat')/sqrt(N)*xhat_dy;
    NMSE_iter_final(Iter_in_max,Iter_ref_max) = norm(z0-zhat_dy)/norm(z0);

    if debug
        dbg_opt.threshold = threshold;
        dbg_opt.sigma_w = sigma_w0;
        dbg_opt.xhat = xhat_dy;
    end
    
    %% export variable
    out.NMSE_iter_final = NMSE_iter_final;
    out.omega_hat = omega_hat;
    out.omega_amp = xhat_dy;
    if B==1
        out.omega_amp = xhat_dy*c;
    end
    
end
