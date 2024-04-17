function [F] = F_prim_far(chi,alpha)
% Set initial value to zero
F = 0;
N = length(alpha);
% Get the sum for F
for jj = 1:N
    F = F - jj*alpha(jj)*chi^(-jj-1);
end
end

