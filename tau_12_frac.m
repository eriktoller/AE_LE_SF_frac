function [tau_12] = tau_12_frac(z,z1,z2,L,m,beta,mu,g,rho,rho_w,k,Omega_w_func)
%Tau 11 for uniform seepage flow in a half-space
%   Detailed explanation goes here

% Getting the chi- and Z-coordinates
chi = chi_from_z(z,z1,z2);
chi_bar = conj(chi);
z0 = 0.5*(z1+z2);

% Calculating the series
dphi = 0;
dphi_bar = 0;
beta_bar = conj(beta);
for n = 1:m
    dphi = dphi + beta_bar(n)*n* (chi^(1-n))/(chi^2-1);
    dphi_bar = dphi_bar + beta(n)*n* (chi_bar^(1-n))/(chi_bar^2-1);
end

% Multiplying the constants
dphi = dphi*(4/L)*exp(-1i*mu);
dphi_bar = dphi_bar*(4/L)*exp(1i*mu);

% Getting the body force
Om = Omega(z,g,rho,rho_w,k,Omega_w_func);

% Calcualting the taus
tau_12 = - exp(1i*mu)*dphi - exp(-1i*mu)*dphi_bar + 0*(conj(Om) + Om);

end

