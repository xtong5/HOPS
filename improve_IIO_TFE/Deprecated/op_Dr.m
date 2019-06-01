function [u]=op_Dr(u_n,D,N_theta,d,R)
% D_r = r partial_r
u = R*op_Partial_r(u_n,D,N_theta,d);
end