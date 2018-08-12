function [u]=op_D_theta(u_n,N_r,N_theta,p)
% D_theta = partial_theta
u=zeros(N_theta*(N_r+1),1);
u_p=zeros(N_theta,1);
for i = 1:N_r+1
 for j = 1:N_theta
    u_p(j) = u_n(i+(j-1)*(N_r+1));
 end
 u_theta = ifft((1i*p).*fft(u_p));
 for j = 1:N_theta
    u(i+(j-1)*(N_r+1)) = u_theta(j);
 end 
end