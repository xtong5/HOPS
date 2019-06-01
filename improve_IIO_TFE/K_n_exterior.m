function K_n=K_n_exterior(u_n,N_r,N_theta,f,T_p,d)
u = zeros(N_theta,1);
for j = 1:N_theta
        u(j) = u_n(1+(j-1)*(N_r+1));
end
K_n = fft(f.* ifft(T_p.*fft(u))/d); 
end