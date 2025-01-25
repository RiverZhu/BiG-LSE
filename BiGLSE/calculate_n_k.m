function [eta,kappa] = calculate_n_k(q,delta_epsilon,nu_q)
% s=1
% Calculate p_q
upper_1 = (q+delta_epsilon)./sqrt(nu_q);
lower_1 = (q-delta_epsilon)./sqrt(nu_q);
Prob = 0.5 * erfc( - (1/sqrt(2))*upper_1 ) - 0.5 * erfc( - (1/sqrt(2))*lower_1 ); 
%  Prob2 =  0.5 * erfc((1/sqrt(2))*lower_1 )-0.5 * erfc((1/sqrt(2))*upper_1 );


p_q = 1/(2*delta_epsilon).*Prob;

% Calculate eta
upper_2 = -(delta_epsilon+q).^2./(2*nu_q);
lower_2 = -(delta_epsilon-q).^2./(2*nu_q);
beta = nu_q./(2*sqrt(2*pi*nu_q).*delta_epsilon);
e_diff = exp(upper_2)-exp(lower_2);
combination_1 = beta.*e_diff./p_q;
eta = q+combination_1;

% Calculate kappa

kappa_part1 = q.^2+2*q.*combination_1;
e_add = exp(lower_2).*(-lower_1)+exp(upper_2).*upper_1;
kappa_part2 = -nu_q./(2*sqrt(2*pi)*delta_epsilon).*e_add./p_q+nu_q;
dbg1 = e_add./(2*sqrt(2*pi)*delta_epsilon.*p_q);
dbg2 = nu_q./(2*sqrt(2*pi)*delta_epsilon).*e_add./p_q;
dbg3 = -dbg2 + nu_q; % should be zero but not
% kappa = kappa_part1+kappa_part2;

kappa = kappa_part2-(combination_1).^2;

I = find(Prob <= 1e-8); % deal with the case with extremely small Prob
if ~isempty(I)
    sigma_d = delta_epsilon^2/3;
    eta(I) = nu_q(I).*sigma_d./(sigma_d+nu_q(I));
    kappa(I) = nu_q(I).*(q(I)./nu_q(I));
end


end

