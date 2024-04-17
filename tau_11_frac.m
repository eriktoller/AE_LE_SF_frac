function [tau_11] = tau_11_frac(z,z1,z2,L,m,beta,mu)
%Tau 11 for uniform seepage flow in a half-space
%   Detailed explanation goes here

% Getting the chi- and Z-coordinates
chi = chi_from_z(z,z1,z2);
chi_bar = conj(chi);
Z = (2*z - z1 - z2)/(z2 - z1);
Z_bar = conj(Z);

% Calculating the series
dphi = 0;
ddphi = 0;
dpsi = 0;
beta_bar = conj(beta);
for n = 1:m
    dphi = dphi + beta_bar(n)*n* (chi^(1-n))/(chi^2-1);
    ddphi = ddphi - beta_bar(n)*n*( chi^(2-n)/((chi^2-1)^3)*((n+1)*chi^2 - n + 1) );
    dpsi = dpsi - beta(n)*n* (chi^(1-n))/(chi^2-1);
end

% Multiplying the constants
dphi = dphi*(4/L)*exp(-1i*mu);
ddphi = ddphi*(16/L^2)*exp(-2i*mu);
dpsi = dpsi*(4/L)*exp(-1i*mu);

% Calcualting the taus
tau_11 = - 0.5*L*(Z-Z_bar)*ddphi - exp(-1i*mu)*(dphi+dpsi);

end

