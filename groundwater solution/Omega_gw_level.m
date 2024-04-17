function [Omega] = Omega_gw_level(z,z_ref,z1,z2,a,Q)
chi = chi_from_z(z,z1,z2);
chi_far = chi_from_z(z_ref,z1,z2);
if chi*conj(chi) < (1-1e-6)
    Omega = complex(NaN,NaN);
else
    Omega = 0;
    for nn = 1:length(a)-1
        Omega = Omega + a(nn+1)/chi^(nn);
    end
    
    % Add the discharge
    Omega = Omega + Q/(2*pi)*log(chi/abs(chi_far));
end
end

