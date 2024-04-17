function [T] = T_total(z,z1,z2,L,m,beta,mu,rt,zt,a,b,ang,sigma_11,g,rho,rho_w,k,Omega_w_func,m_not_f,m_not_t)
% Calcuating the contribution from ech element excluding element m_not
[tau_11,tau_12] = tau_total(z,z1,z2,L,m,beta,mu,rt,zt,a,b,sigma_11,g,rho,rho_w,k,Omega_w_func,m_not_f,m_not_t);

T = -0.5i*(tau_11*exp(2i*ang) - tau_12);

end

