function [beta] = AE_crack_solver(ts,tn,term,A)
% Building the matrices
T_s = ts.*term;
T_n = tn.*term;

% Solvning the matrices
b1 = A\T_s;
b2 = A\T_n;

% Solving for alpha and beta
beta = complex(b2,b1); 

end

