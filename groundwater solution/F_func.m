function [F] = F_func(chi,beta,n)
% Get the value for p
p = p_func(chi,beta,n);
% Sum the values for F
F = (chi^n + chi^(-n))*log( (chi-1)/(chi+1) ) + p;
end

