function [beta] = get_beta(m)
% Creat vector with 1/j for j in (1:m)
B = 1./(1:m);
% Get the odd positions = 1 and even = 0
odd_num = mod(1:m,2);
% Mutiply to set all even values equal to zero
beta = B.*odd_num;
end

