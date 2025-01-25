function [lambda,m_t,V_t] = denoise_input_EM(pi,mu,sigvar,rhat,vr)
% Non-linear variable estimation 
M = log(vr./(vr+sigvar))+(abs(rhat)).^2./vr-(abs(rhat-mu)).^2./(vr+sigvar);
lambda = pi./(pi+(1-pi).*exp(-M));
m_t = (rhat.*sigvar+vr.*mu)./(vr+sigvar);
V_t = vr.*sigvar./(vr+sigvar);
end

