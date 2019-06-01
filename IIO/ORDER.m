
clear all;
close all;
warning off;

L = 2*pi;
lambda = 0.45;
n_u = 1;
% n_w = 2.5;
k_0 = L/lambda;
k_u = 1.1;
a= 1;
% k_w = n_w*k_0;
k_w = 2;
eta = 3.4;
Mode =1;
if(Mode==1)
  sigma_u = 1;
  sigma_w = 1;
else
  sigma_u = (0.4*k_u/L)^2;
  sigma_w = (0.4*k_w/L)^2;
end

pp = [1:100]';

Jp_p = diff_besselj(pp,1,k_w*a);
Jp = besselj(pp,k_w*a);
Hp_p = diff_besselh(pp,1,k_u*a);
Hp = besselh(pp,k_u*a);

figure(1)

Y_p = abs(Jp_p./Jp);
A = polyfit(log(pp),log(Y_p),1); slope1 = A(1);
subplot(2,2,1);
loglog(pp,Y_p,'b-o')
title(['Yp: slope=',num2str(slope1)])

Z_p = abs(Hp_p./Hp);
B = polyfit(log(pp),log(Z_p),1); slope2 = B(1);
subplot(2,2,2);
loglog(pp,Z_p,'b-o')
title(['Zp; slope=',num2str(slope2)])

Q_p =abs ((-k_u*sigma_u*a*Hp_p -1i*eta*Hp)./(-k_u*sigma_u*a*Hp_p +1i*eta*Hp));
C = polyfit(log(pp),log(Q_p),1); slope3 = C(1);
subplot(2,2,3);
loglog(pp,Q_p,'b-o');
title(['Q: slope=',num2str(slope3)])

S_p = abs((k_w*sigma_w*a*Jp_p -1i*eta*Jp)./(k_w*sigma_w*a*Jp_p +1i*eta*Jp));
D = polyfit(log(pp),log(S_p),1); slope4 = D(1);
subplot(2,2,4);
loglog(pp,S_p,'b-o');
title(['S: slope=',num2str(slope4)])

figure(2)
Delta = abs( -k_w*a*Hp.*Jp_p + k_u*a*Hp_p.*Jp );
F = polyfit(log(pp),log(Delta),1); slope5 = F(1);
semilogy(pp,Delta,'b-o');
title(['slope=', num2str(slope5)])
