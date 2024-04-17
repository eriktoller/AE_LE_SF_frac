function [w] = w_uni(z,sigma_11,kappa,G)
% Getting the z- and Z-coordinates
z_bar = conj(z);

% Calculating the series
phi_bar = -0.5*sigma_11*z_bar;
dphi = -0.5*sigma_11;
psi = -0.5*sigma_11*z;

% Calculating w
w =  (z-z_bar)*dphi + kappa*phi_bar + psi ;
end

