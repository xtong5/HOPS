clear all
close all

L = 2*pi;
M = 11;
lambda = linspace(0.3,0.8,M);
% lambda = 0.6;
k_zero = 2*pi./lambda;
n_u = 2.5;
m_w = zeros(M,1);
epsilon_w = zeros(M,1);
epsilon_u_plot = zeros(M,1);
epsilon_w_plot = zeros(M,1);
epsilon_u = n_u^2;
for i = 1:M
    [n_w(i),epsilon_w(i)] = ri_perm(lambda(i),'SILVER');
    epsilon_u_plot(i) = epsilon_u;
    epsilon_w_plot(i) = epsilon_w(i);
end

[min,r] = min(abs(epsilon_u - real(-epsilon_w)));
lambda_crit = lambda(r);

k_u = n_u*k_zero; 
k_w = n_w.*k_zero; 
RunNumber = 0;
Mode = 2; 
% Mode = 1;

if(RunNumber==0)
  N_theta = 16;
%   N_theta = 16;
  a = 0.025;
  b = 10*a;
  N = 16;
  Eps = 0.01*a;
elseif(RunNumber==1)
  % Small Deformation
  Eps = 0.02;
  N_theta = 64;
  a = 1; 
  b = 5;
  N = 16;
elseif(RunNumber==2)
  % Big Deformation (inside disk)
  Eps = 0.3;
  N_theta = 64;
  a = 2.0;
  b = 5;
  N = 16;
elseif(RunNumber==3)
  % Big Deformation (outside disk)
  Eps = 0.75;
  N_theta = 64;
  a = 2.0;
  b = 5;
  N = 16;
end

fprintf('test\n');
fprintf('-------------\n');
fprintf('RunNumber = %d\n',RunNumber);
fprintf('Eps = %g  a = %g  b = %g\n',Eps,a,b);
fprintf('N_theta = %d N = %d M = %d\n',N_theta,N,M);
fprintf('\n');

theta = (L/N_theta)*[0:N_theta-1]';
p = [0:N_theta/2-1,-N_theta/2:-1]';

f = exp(cos(theta));
f_theta = ifft(1i*p.*fft(f));
A = a+Eps*f;

P_0 = 1;

BU = zeros(N_theta,M);
U = zeros(N_theta,M);
W = zeros(N_theta,M);
B_far = zeros(N_theta,N+1);
BU_norm = zeros(M,1);
U_norm = zeros(M,1);
W_norm = zeros(M,1);
U_n = zeros(N_theta,N+1,M);
W_n = zeros(N_theta,N+1,M); 

for m =1:M

zeta = -P_0*exp(-1i*k_u(m).*A.*sin(theta));
psi = (1i*k_u(m)).*A.*zeta;
zeta_n = zeros(N_theta,N+1);
psi_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);
Sin = sin(theta);
Exp = exp(-1i*k_u(m)*a.*Sin);
zeta_n(:,0+1) = -Exp; 
psi_n(:,0+1) = (1i*k_u(m))*a*Sin.*Exp;
for n=1:N
f_n = f.*f_n/n;
if n > 1
   f_nmo = f.*f_nmo/(n-1);
end

zeta_n(:,n+1) = -Exp.*(-1i*k_u(m))^n.*f_n.*Sin.^n;
psi_n(:,n+1) = (1i*k_u(m)).*Exp.*(a*(-1i*k_u(m))^n.*f_n.*Sin.^(n+1)+...
    (f.*Sin-f_theta.*cos(theta)).*(-1i*k_u(m))^(n-1).*f_nmo.*Sin.^(n-1));
end

if(Mode==1)
  tau2 = 1;
else
  tau2 = k_u(m)^2/k_w(m)^2;
end


   U_n(:,:,m) = twolayer_dno_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,tau2,...
    p,k_u(m),k_w(m),a,N_theta,N);
   apn_fe = field_fe_helmholtz_polar_exterior(U_n(:,:,m),f,k_u(m),a,p,N_theta,N);
   W_n(:,:,m) = U_n(:,:,m) - zeta_n;
   
   for n=1:N+1
     B_far(:,n)=ifft(apn_fe(:,n).*besselh(p,k_u(m).*b)./besselh(p,k_u(m).*a));
   end 

   for j=1:N_theta
     BU(j,m) = taylorsum(B_far(j,:),Eps,N);
     U(j,m) = taylorsum(U_n(j,:,m),Eps,N);
     W(j,m) = taylorsum(W_n(j,:,m),Eps,N);
   end
   
   
   
   U_norm(m) = norm(U(:,m),2)/sqrt(N_theta);
   BU_norm(m) = norm(BU(:,m),2)/sqrt(N_theta);
   W_norm(m) = norm(W(:,m),2)/sqrt(N_theta);
   
end
% plot(lambda,U_norm,'b-*',lambda,BU_norm,'r-o')
% plot(lambda,error_BU,'b-*',lambda,error_U,'r-o')
% plot(theta,real(BU),'b-*',theta,real(BU_exact),'r-o')
% plot(theta,U,'b-*',theta,U_exact,'r-o')
% plot(theta,error_U,'b-*',theta,error_BU,'r-o')
% plot(lambda,U_norm,'b-o',lambda,W_norm,'g-*',lambda,BU_norm,'c-d')



lambda_crit_plot = lambda_crit*ones(M,1);
norm_max = max([max(U_norm),max(W_norm),max(BU_norm)]);
yy = linspace(0,norm_max,M)';

figure(1);
plot(lambda,U_norm,'b-o',lambda,W_norm,'g-*',...
    lambda,BU_norm,'c-d',lambda_crit_plot,yy,'r--');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|U|_2$ and $|W|_2$','interpreter','latex');
title('$|U|_2$ and $|W|_2$ versus $\lambda$','interpreter','latex');
legend('|U|_2','|W|_2','|BU|_2','lambda_c');



%
% Plot of critical lambda
%

eps_max = max(max(epsilon_u_plot),max(-real(epsilon_w_plot)));
yy = linspace(0,eps_max,M)';

figure(2);
plot(lambda,epsilon_u_plot,'b-o',lambda,-real(epsilon_w_plot),'g-*',...
    lambda_crit_plot,yy,'r--');
xlabel('$\lambda$','interpreter','latex');
ylabel('$\epsilon$','interpreter','latex');
title('$\epsilon_u$ and -Real($\epsilon_w$) versus $\lambda$','interpreter','latex');
legend('epsilon_u','-Re[epsilon_w]','lambda_c');

% save('Perturbation.mat','M','lambda','U_norm','W_norm','BU_norm',...
%     'lambda_crit','epsilon_u_plot','epsilon_w_plot')




 

