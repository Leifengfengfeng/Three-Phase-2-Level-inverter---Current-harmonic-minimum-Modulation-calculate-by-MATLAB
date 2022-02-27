function [c,ceq] = ncf_B(alpha)
global M;

c = [];
ceq = Um_B(alpha,1) - M; 

end
