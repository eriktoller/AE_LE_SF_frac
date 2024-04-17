function [W] = W_tunnel(z,rt,zt,t,Q)
Z = bigz_from_z(z,rt,zt);
if Z*conj(Z) < (1-1e-6)
    W = complex(NaN,NaN);
else
    W = 0;
    for nn = 1:length(t)-1
        W = W + nn*t(nn+1)/Z^(nn+1);
    end
    
    % Add the discharge
    W = W - Q/(2*pi)*(1/Z);
end
W = W/rt;
end

