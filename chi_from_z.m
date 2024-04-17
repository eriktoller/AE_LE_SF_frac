function [chi] = chi_from_z(z,z1,z2)
% Mapp small z to big Z
Z = (2*z - z1 - z2)/(z2 - z1);
% Mapp Z to chi
chi = Z + sqrt(Z-1)*sqrt(Z+1);
end

