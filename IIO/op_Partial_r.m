function [u]=op_Partial_r(u_n,D,N_theta,d)
% partial_r 
D_r = kron(eye(N_theta),D*(2/d));
u = D_r*u_n;
end