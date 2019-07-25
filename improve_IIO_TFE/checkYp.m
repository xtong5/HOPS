L = 2*pi;
lambda = 0.45;
n_u = 1;
n_w = 2.5;
k_0 = L/lambda;
k_u = n_u*k_0; 
k_w = n_w*k_0;

N_theta = 64; 
p = [0:N_theta/2-1,-N_theta/2:-1]';

a = 0.025;
c = 0.1*a;

Y_p = k_w.*diff_besselj(p,1,k_w.*c)./besselj(p,k_w.*c);
Z_p = k_u.*diff_besselh(p,1,k_u.*b)./besselj(p,k_u.*b);