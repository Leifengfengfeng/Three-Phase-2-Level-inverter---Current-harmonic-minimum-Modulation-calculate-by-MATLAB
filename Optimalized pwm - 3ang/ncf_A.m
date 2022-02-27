function [c,ceq] = ncf_A(alpha)
global M;

c = [];
ceq = Um_A(alpha,1) - M; 

end
