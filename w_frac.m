function [w] = w_frac(z,z1,z2,L,m,beta,mu,kappa,G)
% Getting the chi- and Z-coordinates
chi = chi_from_z(z,z1,z2);
chi_bar = conj(chi);
z0 = 0.5*(z1+z2);
Z = (2*z - z1 - z2)/(z2 - z1);
Z_bar = conj(Z);

% Calculating the series
phi_bar = 0;
dphi = 0;
psi = 0;
beta_bar = conj(beta);
for n = 1:m
    dphi = dphi + beta_bar(n)*n*(4/L)*exp(-1i*mu) * (chi^(1-n))/(chi^2-1);
    phi_bar = phi_bar - beta(n)*chi_bar^(-n);
    psi = psi + beta(n)*chi^(-n);
end

% Calculating w
w =  0.5*L*(Z-Z_bar)*dphi + exp(-1i*mu)*kappa*phi_bar + exp(-1i*mu)*psi ;

end

