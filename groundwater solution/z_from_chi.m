function [z] = z_from_chi(chi,z1,z2)
% Mapp Z from chi
Z = .5.*(chi+1./chi);
% Mapp small z from big Z
z = .5.*( Z.*(z2-z1)+z1+z2 );
end

