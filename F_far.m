function [F] = F_far(chi,alpha)
% Set the initial value to zero
F = 0;
N = length(alpha);
% Get the sum for F
for jj = 1:N
    F = F + alpha(jj)*chi^(-jj);
end
end

