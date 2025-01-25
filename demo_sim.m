%% Initialization
clc  
clear
close all
addpath('.\fuc')
addpath('.\superfast-LSE')
addpath('.\GNOMP_K')
addpath('.\BiGLSE')
seed = randi(1000);
rng(seed)

N = 1024; % number of measurement
B = 2; % quantization bit, could be selected in {1,2,3,\cdots, inf}
idx = 1; % in {1,2,3,4} to choose different scenario
% generate line spectral signal
switch idx
    case 1 % frequency interval is specified, little number of targets scenario
        d = 2*pi/N;
        parameter.K_MAX = 30; K = 10;
    case 2 % frequency interval is not limited, little number of targets scenario
        d = 0;
        parameter.K_MAX = 30; K = 10;
    case 3 % frequency interval is specified, multi targets scenario
        d = 2*pi/N;
        parameter.K_MAX = 200; K = 100;
    case 4 % frequency interval is not limited, multi targets scenario
        d = 0;
        parameter.K_MAX = 200; K = 100;
end

% line spectral signal parameter initialize
pivar = K/N; mu_x = 0; sigvar = 1e0; L = 4*N;
SNR_min = 25; SNR_max = 25; % SNR 
SNR = SNR_min*ones(K,1) + rand(K,1)*(SNR_max - SNR_min);
SNR(1) = 20;
theta = zeros(K,1);
theta(1) = pi*2*rand;
for k = 2:K
    th = pi * 2*rand;
    while min(abs((wrapToPi(th-theta(1:k-1))))) < d
        th = pi * 2*rand;
    end
    theta(k) = th;
end
A = exp(1j*(0:N-1).'*theta.')/sqrt(N);
sigma_w = 1;
noise = sqrt(sigma_w/2)*(randn(N,1)+1j*randn(N,1));
r = sqrt(10.^(SNR/10)*sigma_w);
x = r.*exp(1j*2*pi*rand(K,1));
z0 = A*x;

if B == inf
    y = z0 + noise;
else
    y_unq = z0 + noise;
    y = [real(y_unq);imag(y_unq)];
    nbins = 2^B;
    sig_z = L*pivar*(abs(mu_x)^2 + sigvar)/N;
    y_max = 3*sqrt(sig_z); y_min = -y_max;
    alpha = (y_max - y_min)/(nbins);
    yy = floor((y-y_min)/alpha);
    index1 = find(y>=y_max); yy(index1) = nbins-1;
    index2 = find(y<y_min); yy(index2) = 0;
    y = yy;
end

%% algorithm run
if B==inf
    % (1) SuperFast-LSE
    t_in = tic;
    out_sf = superfast_lse(y, (1:N)', N, 'verbose',false, 'plot',false);
    time_sf = toc(t_in);
    omega_hat_sf = wrapTo2Pi(-out_sf.tau*2*pi);
    [omega_hat_sf, od] = sort(omega_hat_sf);
    omega_amp_sf = out_sf.alpha(od)*sqrt(N);
    NMSE_sf = 20*log10(norm(z0-out_sf.h)/norm(z0));
    
    % (2) NOMP
    p_fa = 1e-2;
    tau = sigma_w * ( log(N) - log( log(1/(1-p_fa)) ) );
    S = eye(N);
    t_in = tic;
    [omega_hat_nomp, omega_amp_nomp, residue] = extractSpectrum(y, S, tau);
    time_nomp = toc(t_in);
    zhat_nomp = exp(1j*((0:N-1)')*omega_hat_nomp')/sqrt(N)*omega_amp_nomp;
    NMSE_nomp = 20*log10(norm(z0-zhat_nomp)/norm(z0));
    
    % caculate algorithm lowerbound
    SNR_ndB = norm(z0)^2/(N*sigma_w);
    tmp = A*pinv(A);
    NMSE_ora = 10*log10(trace(tmp'*tmp)/(N*SNR_ndB));
    
end

% (3) BiGLSE(para.)
OS = 3; L_prl = OS*N;
if B==inf
    Init_prl = LSEOpt('N_len',N, 'L_len',L_prl, 'Bit',B, 'ref_loop',1, 'iter_max',50, ...
        'Iter_ref_max',6, 'Iter_in_max',50, 'maxBadSteps',5);
else
    Init_prl = LSEOpt('N_len',N, 'L_len',L_prl, 'Bit',B, 'ref_loop',1, 'iter_max',50, 'Iter_ref_max',6,...
        'Iter_in_max',50, 'maxBadSteps',5,'y_min',y_min, 'alpha',alpha, 'sigma_w',sigma_w);
end
theta_g = (0:2*pi/L_prl:2*pi/L_prl*(L_prl-1)); a_theta = exp(1j*(0:N-1).'*theta_g)/sqrt(N);
t_in = tic;
if B == inf
    [out_bl_prl, dbg_opt_prl] = BiGLSE_prl(y, z0, a_theta, Init_prl, []);
else
    [out_bl_prl, dbg_opt_prl] = BiGLSE_prl_lq(y, z0, a_theta, Init_prl, []);
end
time_bl_prl = toc(t_in);
NMSE_iter_final_prl = 20*log10(out_bl_prl.NMSE_iter_final);
NMSE_prl = NMSE_iter_final_prl(end);
omega_hat_prl = out_bl_prl.omega_hat; omega_amp_prl = (out_bl_prl.omega_amp).*exp(-1j*(N-1)/2);

% (4) BiGLSE(seri.)
OS = 4; L_seq = OS*N;
parameter.P_fa = 1e-6; parameter.OverSampling = OS;
if B==inf
    Init_seq = LSEOpt('N_len',N, 'L_len',L_seq, 'Bit',B, 'Iter_ref_max',7, 'Iter_in_max',75);
else
    Init_seq = LSEOpt('N_len',N, 'L_len',L_seq, 'Bit',B, 'Iter_ref_max',7, 'Iter_in_max',75,...
        'y_min',y_min, 'alpha',alpha, 'sigma_w',sigma_w);
end
t_in = tic;
if B == inf
    [out_bl_seq, dbg_opt_seq] = BiGLSE_seq(y, z0, parameter, Init_seq, []);
else
    [out_bl_seq, dbg_opt_seq] = BiGLSE_seq_lq(y, z0, parameter, Init_seq, []);
end
time_bl_seq = toc(t_in);
NMSE_iter_final_seq = 20*log10(out_bl_seq.NMSE_iter_final);
NMSE_seq = NMSE_iter_final_seq(end);
omega_hat_seq = out_bl_seq.omega_hat; omega_amp_seq = out_bl_seq.omega_amp;

% (5) BiGLSE(vm_para.)
OS = 3; L_vm2 = OS*N;
if B==inf
    Init_vm2 = LSEOpt('N_len',N, 'L_len',L_vm2, 'Bit',B, 'iter_max',50, ...
        'Iter_ref_max',1, 'Iter_in_max',100, 'maxBadSteps',5, 'sigma_w',sigma_w);
else
    Init_vm2 = LSEOpt('N_len',N, 'L_len',L_vm2, 'Bit',B, 'iter_max',50, ...
        'Iter_ref_max',1, 'Iter_in_max',100, 'maxBadSteps',5, ...
        'y_min',y_min, 'alpha',alpha, 'sigma_w',sigma_w);
end
theta_g = 0:2*pi/L_vm2:2*pi/L_vm2*(L_vm2-1); a_theta = exp(1j*(0:N-1).'*theta_g)/sqrt(N);
t_in = tic;
[out_vm2, ~] = BiGLSE_vm_prl(y, z0, a_theta, Init_vm2, []);
time_vm2 = toc(t_in);
NMSE_iter_final_vm2 = 20*log10(out_vm2.NMSE_iter_final);
NMSE_vm2 = NMSE_iter_final_vm2(end);
omega_hat_vm2 = out_vm2.omega_hat;
omega_amp_vm2 =  out_vm2.omega_amp;
% end

% (6) BiGLSE(vm_seri.)
OS = 4; L_vm = OS*N;
parameter.P_fa = 1e-6; parameter.OverSampling = OS;
if B == inf
    Init_vm = LSEOpt('N_len',N, 'L_len',L_vm, 'Bit',B, 'Iter_ref_max',1, 'Iter_in_max',100, ....
        'sigma_w', sigma_w);
else
    Init_vm = LSEOpt('N_len',N, 'L_len',L_vm, 'Bit',B, 'Iter_ref_max',1, 'Iter_in_max',100, ....
        'y_min',y_min, 'alpha',alpha, 'sigma_w',sigma_w);
end
dbg_ipt.theta = theta; dbg_ipt.x = x;
t_in = tic;
[out_vm, dbg_opt_vm] = BiGLSE_vm_seq(y, z0, parameter, Init_vm, dbg_ipt);
time_vm = toc(t_in);
NMSE_iter_final_vm = 20*log10(out_vm.NMSE_iter_final);
NMSE_vm = NMSE_iter_final_vm(end);
omega_hat_vm = out_vm.omega_hat; omega_amp_vm = out_vm.omega_amp;

%% plot result
figure(1)
plot(1:Init_vm.Iter_in_max, NMSE_iter_final_vm, 'r');
hold on
plot(1:Init_vm2.Iter_in_max, NMSE_iter_final_vm2, 'c');
plot(1:Init_seq.Iter_in_max, NMSE_iter_final_seq(:,end), 'k');
plot(1:Init_prl.Iter_in_max, NMSE_iter_final_prl(:,end), 'm');
legend('BiG-LSE(vm\_seri.)', 'BiG-LSE(vm\_para.)', 'BiG-LSE(seri.)', 'BiG-LSE(para.)')
xlabel('number of iterations')
ylabel('NMSE (dB)')


figure(2)
stem(theta, abs(x), 'bo')
hold on
if B == inf
    stem(omega_hat_sf, abs(omega_amp_sf), 'rx')
    stem(omega_hat_nomp, abs(omega_amp_nomp), 'c>')
end
stem(omega_hat_prl, abs(omega_amp_prl), 'm+')
stem(omega_hat_seq, abs(omega_amp_seq), 'ks')
if B ==inf
    legend('True', 'Superfast LSE', 'NOMP', 'BiG-LSE(para.)', 'BiG-LSE(seri.)')
else
    legend('True', 'BiG-LSE(para.)', 'BiG-LSE(seri.)')
end

figure(3) 
stem(theta, abs(x), 'bo')
hold on
stem(omega_hat_vm, abs(omega_amp_vm), 'rx')
stem(omega_hat_vm2, abs(omega_amp_vm2), 'c>')
stem(omega_hat_seq, abs(omega_amp_seq), 'ks')
stem(omega_hat_prl, abs(omega_amp_prl), 'm+')
legend('True', 'BiG-LSE(vm\_seri.)', 'BiG-LSE(vm\_para.)', 'BiG-LSE(seri.)', 'BiG-LSE(para.)')
