function [u]=op_A(u_n,f,r,b)
% A = A(f)
A = diag(kron(f,b-r));
u = A*u_n;
end