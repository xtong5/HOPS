% XT 4/19

clear all;
close all;
warning off;

Mode = 2; 

M = 1000;
N_theta = 64;
a = 0.025;

L = 2*pi;
lambda = linspace(0.3,0.8,M);
k_zero = (2*pi./lambda).';
OUT = 'VACUUM';
IN = 'SILVER';

theta = (L/N_theta)*[0:N_theta-1]';
p = [0:N_theta/2-1,-N_theta/2:-1]';

n_u = zeros(M,1);
n_w = zeros(M,1);
epsilon_u = zeros(M,1);
epsilon_w = zeros(M,1);
epsilon_u_plot = zeros(M,1);
epsilon_w_plot = zeros(M,1);

for i = 1:M
    [n_u(i),epsilon_u(i)] = ri_perm(lambda(i),OUT);
    [n_w(i),epsilon_w(i)] = ri_perm(lambda(i),IN);
    epsilon_u_plot(i) = epsilon_u(i);
    epsilon_w_plot(i) = epsilon_w(i);
end

[mini,r] = min(abs(epsilon_u - real(-epsilon_w)));
lambda_crit = lambda(r);

k_u = n_u.*k_zero; 
k_w = n_w.*k_zero; 


if(Mode==1)
  tau2 = ones(M,1);
else
  tau2 = k_u.^2/k_w.^2;
end

pp = 1;

Hp = besselh(pp,k_u.*a);
Jp = besselj(pp,k_w.*a);
Hp_prime = diff_besselh(pp,1,k_u.*a);
Jp_prime = diff_besselj(pp,1,k_w.*a);
Delta_p = -tau2*(k_w.*a).*Hp.*Jp_prime + (k_u.*a).*Hp_prime.*Jp;

lambda_crit_plot = lambda_crit*ones(M,1);
yy = linspace(min(real(Delta_p)),max(real(Delta_p)),M)';
plot(lambda,real(Delta_p))
hold on
plot(lambda_crit_plot,yy,'r--');