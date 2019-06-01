function u_n=FT_inv(u_n_hat,N_r,N_theta)
f_n = zeros(N_theta,1);
u_n = zeros(size(u_n_hat));
for i = 1:N_r+1
    for j = 1:N_theta
        f_n(j) = u_n_hat(i+(j-1)*(N_r+1));
    end
    f_n = ifft(f_n);
    for j = 1:N_theta
        u_n(i+(j-1)*(N_r+1)) = f_n(j);
    end 
end
end