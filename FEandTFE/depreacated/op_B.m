function [u]=op_B(u_n,f_theta,r,b)
% B = B(f)
B = diag(kron(f_theta,b-r));
u = B*u_n;
end