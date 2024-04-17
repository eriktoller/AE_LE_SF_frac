function [p] = p_prim_func(chi,beta,n)
% Set the initial value to zero
p = 0;
if n > 0
    % Cacluate the sum for p'
    for jj = 1:(n-1)
        p = p + jj*beta(n-jj)*(chi^(jj-1) - chi^(-jj-1));
    end
    p = p*2;
end
end

