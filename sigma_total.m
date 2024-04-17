function [sigma_11,sigma_22,sigma_12] = sigma_total(z,z1,z2,L,m,beta,mu,sigma_11,g,rho,rho_w,k,Omega_w_func,m_not)

% Calcuating the contribution from ech element excluding element m_not
[tau_11,tau_12] = tau_total(z,z1,z2,L,m,beta,mu,sigma_11,g,rho,rho_w,k,Omega_w_func,m_not);

% Assigning the sigma
S1 = .5*(tau_11+tau_12);
S2 = .5*(-tau_11+tau_12);
sigma_11 = real(S1);
sigma_22 = real(S2);
sigma_12 = -imag(S1);
end

