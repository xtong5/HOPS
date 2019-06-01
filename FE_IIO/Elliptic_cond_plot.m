% XT 4/19

clear all;
close all;
warning off;

Mode = 2; 

N_theta = 64;
a = 0.025;
b = 10*a;
c = a/10;

L = 2*pi;
% lambda = linspace(0.3,0.8,M);
lambda = 0.45;
% k_zero = (2*pi./lambda).';
k_zero = L/lambda;
OUT = 'VACUUM';
IN = 'SILVER';

theta = (L/N_theta)*[0:N_theta-1]';
p = [0:N_theta/2-1,-N_theta/2:-1]';

% n_u = zeros(M,1);
% n_w = zeros(M,1);
% epsilon_u = zeros(M,1);
% epsilon_w = zeros(M,1);
% epsilon_u_plot = zeros(M,1);
% epsilon_w_plot = zeros(M,1);

% for i = 1:M
%     [n_u(i),epsilon_u(i)] = ri_perm(lambda(i),OUT);
%     [n_w(i),epsilon_w(i)] = ri_perm(lambda(i),IN);
%     epsilon_u_plot(i) = epsilon_u(i);
%     epsilon_w_plot(i) = epsilon_w(i);
% end

[n_u,epsilon_u] = ri_perm(lambda,OUT);
[n_w,epsilon_w] = ri_perm(lambda,IN);

% [mini,r] = min(abs(epsilon_u - real(-epsilon_w)));
% lambda_crit = lambda(r);

k_u = n_u.*k_zero; 
k_w = n_w.*k_zero; 


if(Mode==1)
  sigma_u = 1;
  sigma_w = 1;
else
  sigma_u = 1/epsilon_u;
  sigma_w = 1/epsilon_w;
end

Hp = besselh(p,k_u.*b);
Jp = besselj(p,k_w.*c);
Hp_prime = diff_besselh(p,1,k_u.*b);
Jp_prime = diff_besselj(p,1,k_w.*c);

B_real = real(k_u.*Hp_prime./Hp); B_im = imag(k_u.*Hp_prime./Hp);
A_real = real(k_w.*Jp_prime./Jp); A_im = imag(k_w.*Jp_prime./Jp);

AA = all(A_real>0); BB = all(B_real<0); 
fprintf("real(k_w)= %g \n",real(k_w));
fprintf("Every component of A is nonnegative: %d \n",AA);
fprintf("Every component of B is nonpositive: %d \n",BB);


% plot(p,(p+1)/b,'r-',p,B,'b-')

