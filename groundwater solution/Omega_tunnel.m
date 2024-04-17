function [Omega] = Omega_tunnel(z,z_ref,rt,zt,t,Q)
Z = bigz_from_z(z,rt,zt);
Z_far = bigz_from_z(z_ref,rt,zt);
if Z*conj(Z) < (1-1e-6)
    Omega = complex(NaN,NaN);
else
    Omega = 0;
    for nn = 1:length(t)-1
        Omega = Omega + t(nn+1)/Z^(nn);
    end
    
    % Add the discharge
    Omega = Omega + Q/(2*pi)*log(Z/abs(Z_far));
end
end

