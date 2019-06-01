
clear all;
close all;
warning off;

L = 2*pi;
lambda = 0.45;
n_u = 1;
% n_w = 2.5;
k_0 = L/lambda;
k_u = 1.1;
a=1;
% k_w = n_w*k_0;
k_w = 2;
eta = 3.4;
N_theta = 128;
Mode =1;
if(Mode==1)
  sigma_u = 1;
  sigma_w = 1;
else
  sigma_u = (lambda*k_u/L)^2;
  sigma_w = (lambda*k_w/L)^2;
end

pp = [0:N_theta/2-1,-N_theta/2:-1]';


Jp_p = diff_besselj(pp,1,k_w*a);
Jp = besselj(pp,k_w*a);
Hp_p = diff_besselh(pp,1,k_u*a);
Hp = besselh(pp,k_u*a);

Y_p = abs(Jp_p./Jp);
Z_p = abs(Hp_p./Hp);

S_p = (k_w*sigma_w*a*Jp_p -1i*eta*Jp)./(k_w*sigma_w*a*Jp_p +1i*eta*Jp);
Q_p = (-k_u*sigma_u*a*Hp_p -1i*eta*Hp)./(-k_u*sigma_u*a*Hp_p +1i*eta*Hp);

Root = sqrt(S_p.*Q_p);
Eigen1 = 1+Root;
Eigen2 = 1-Root;

figure(1)
plot(Eigen1,'r*')
hold on
plot(Eigen2,'bo')

figure(2)
plot(S_p,'r*')


figure(3)
plot(Q_p,'r*')
