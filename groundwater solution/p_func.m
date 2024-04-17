function [p] = p_func(chi,beta,n)
% Set the initial value to zero
p = 0;
if n > 0
    p = beta(n);
    % Calculate the sum for p
    for jj = 1:(n-1)
        p = p + beta(n-jj)*(chi^(jj)+chi^(-jj));
    end
    p = p*2;
end
end

