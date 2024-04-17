function [alpha] = get_alpha(beta,n,mfar)
% Set all initial value of alpha to zero
alpha = zeros(1,mfar);
% Assing alpha for j < n-1
for jj = 1:(n-1)
    alpha(jj) = 2*(beta(n-jj) - beta(n+jj));
end
% Assing alpha för j > n+1
for jj = (n+1):mfar
    alpha(jj) = -2*(beta(jj-n) + beta(n+jj));
end
end

