function [tau_11,tau_12] = tau_uni_p(sigma_11)
% Calculating the phi and psi
dphi = -.5*sigma_11;
dpsi = dphi;

% Calcualting the taus
tau_11 = - dphi - dpsi;
tau_12 = - dphi - dphi;
end

