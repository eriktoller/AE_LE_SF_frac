function [W] = W_total(z,W_uni,z1b,z2b,b,z1a,z2a,a,beta,alpha,m,rt,zt,t,Q,z1gw,z2gw,c,m_notb,m_nota,m_nott,m_notgw)
% Get W for uniform flow
W = W_uni;
% Add the W for the draining fractures (which are not number m_not)
for ii = 1:length(z1b)
    if ne(ii,m_notb)
        chi = chi_from_z(z,z1b(ii),z2b(ii));
        W = W + W_dp(chi,z1b(ii),z2b(ii),b(ii,:),beta,alpha,m);
    end
end
% Add the W for the blocking fractures (which are not number m_not)
for ii = 1:length(z1a)
    if ne(ii,m_nota)
        chi = chi_from_z(z,z1a(ii),z2a(ii));
        W = W + W_db(chi,z1a(ii),z2a(ii),a(ii,:),beta,alpha,m);
    end
end
% Add W for the tunnel
for ii = 1:length(zt)
    if ne(ii,m_nott)
        W = W + W_tunnel(z,rt(ii),zt(ii),t(ii,:),Q(ii));
    end
end
% Add W for groudnwater level
for ii = 1:length(z1gw)
    if ne(ii,m_notgw)
        W = W + W_gw_level(z,z1gw(ii),z2gw(ii),c(ii,:),Q(length(zt)+ii));
    end
end
end

