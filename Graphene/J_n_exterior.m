function J_n=J_n_exterior(u_n,N_r,N_theta,f,T_p,d)
u = zeros(N_theta,1);
for j = 1:N_theta
        u(j) = u_n(1+(j-1)*(N_r+1));
end
J_n = d/2*fft(-f.* ifft(T_p.*fft(u))/d); 
end