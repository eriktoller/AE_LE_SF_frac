function [Omega] = Omega_total(z,Phi0,W_uni,z1b,z2b,b,z1a,z2a,a,beta,alpha,m,z_ref,rt,zt,t,Q,z1gw,z2gw,c,m_notb,m_nota,m_nott,m_notgw)
% Get the Omega for uniform flow
Omega = Omega_uni(z,W_uni) + Phi0;
% Add Omega for the draining fractures
for ii = 1:length(z1b)
    if ne(ii,m_notb)
        chi = chi_from_z(z,z1b(ii),z2b(ii));
        Omega = Omega + Omega_dp(chi,b(ii,:),beta,alpha,m);
    end
end
% Add Omega for the blocking fractures
for ii = 1:length(z1a)
    if ne(ii,m_nota)
        chi = chi_from_z(z,z1a(ii),z2a(ii));
        Omega = Omega + Omega_db(chi,a(ii,:),beta,alpha,m);
    end
end
% Add Omega for the tunnel
for ii = 1:length(zt)
    if ne(ii,m_nott)
        Omega = Omega + Omega_tunnel(z,z_ref,rt(ii),zt(ii),t(ii,:),Q(ii));
    end
end
% Add Omega for groudnwater level
for ii = 1:length(z1gw)
    if ne(ii,m_notgw)
        Omega = Omega + Omega_gw_level(z,z_ref,z1gw(ii),z2gw(ii),c(ii,:),Q(length(zt)+ii));
    end
end
end

