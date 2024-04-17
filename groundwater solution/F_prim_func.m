function [F] = F_prim_func(chi,beta,n)
% Get teh value for p'
p = p_prim_func(chi,beta,n);
% Sum teh value for F'
F = n*(chi^(n-1)-chi^(-n-1))*log( (chi-1)/(chi+1) ) + ...
    (chi^n + chi^(-n))*2/(chi^2-1) + p;
end

    