function [sigma_1,sigma_2,theta_p] = major_minor_stress(z,tau_func)
% Get the tau_11 and tau_12
[tau_11,tau_12] = tau_func(z);

% Calculate the sigma
S1 = .5*(tau_11+tau_12);
S2 = .5*(-tau_11+tau_12);
sigma_11 = real(S1);
sigma_22 = real(S2);
sigma_12 = -imag(S1);

% Calculate the terms for the prinicpal stresses
frac1 = (sigma_11 + sigma_22)/2;
sqrt1 = sqrt( ((sigma_11 - sigma_22)/2)^2 + sigma_12^2 );

% Calcualte the major and minor principal stress
sigma_1 = frac1 + sqrt1;
sigma_2 = frac1 - sqrt1;

% Calcualte the angel of the major pricipal stress
theta_p = -0.5*imag(log(tau_11));
end

