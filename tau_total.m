function [tau_11,tau_12] = tau_total(z,z1,z2,L,m,beta,mu,rt,zt,a,b,sigma_11,g,rho,rho_w,k,Omega_w_func,m_not_f,m_not_t)

tau_11f = 0;
tau_12f = 0;
for ii = 1:length(z1)
    if ne(ii,m_not_f)
        tau_11f = tau_11f + tau_11_frac(z,z1(ii),z2(ii),L(ii),m,beta(ii,:),mu(ii));
        tau_12f = tau_12f + tau_12_frac(z,z1(ii),z2(ii),L(ii),m,beta(ii,:),mu(ii),g,rho,rho_w,k,Omega_w_func);
    end
end

tau_11t = 0;
tau_12t = 0;
for ii = 1:length(zt)
    if ne(ii,m_not_t)
        tau_11t = tau_11t + tau_11_tunnel(z,rt,zt,a,b);
        tau_12t = tau_12t + tau_12_tunnel(z,rt,zt,a,b,g,rho,rho_w,k,Omega_w_func);
    end
end

% Add Omega as a separate function
Om = 0;
if m_not_f + m_not_t == 0
    Om = Omega(z,g,rho,rho_w,k,Omega_w_func);
end

[tau_11p,tau_12p] = tau_uni_p(sigma_11);
tau_11 = tau_11f + tau_11t + tau_11p;
tau_12 = tau_12f + tau_12t + tau_12p + conj(Om) + Om;
end

