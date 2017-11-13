function u_n_hat=FT(u_n,N_r,N_theta)
f_hat = zeros(N_theta,1);
u_n_hat = zeros(size(u_n));
for i = 1:N_r+1
    for j = 1:N_theta
        f_hat(j) = u_n(i+(j-1)*(N_r+1));
    end
    f_hat = fft(f_hat);
    for j = 1:N_theta
        u_n_hat(i+(j-1)*(N_r+1)) = f_hat(j);
    end 
end
end