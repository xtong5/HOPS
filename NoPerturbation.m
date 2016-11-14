
clear all
close all

L = 2*pi;
%lambda = [0.2:0.005:1];
M = 401;
lambda = linspace(0.3,0.8,M);
% lambda = 0.6;
k_zero = 2*pi./lambda;
n_u = 2.5;
% n_w = 1.3;
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
Mode = 2; %check 
% Mode = 1;



if(RunNumber==0)
  Eps = 0;
  N_theta = 64;
%   N_theta = 16;
  a = 0.025;
  b = 10*a;
  N = 0;
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

% f = exp(cos(theta));
% f_theta = -sin(theta).*f;
% A = a+Eps*f;

P_0 = 1;


BU_exact = zeros(N_theta,M);
U_exact = zeros(N_theta,M);
W_exact = zeros(N_theta,M);
BU_norm = zeros(M,1);
U_norm = zeros(M,1);
W_norm = zeros(M,1);

for m =1:M

if(Mode==1)
  tau2 = 1;
else
  tau2 = k_u(m)^2/k_w(m)^2;
end

zeta = -P_0*exp(-1i*k_u(m)*a*sin(theta));
psi = -P_0*(-1i*k_u(m)*sin(theta)).*exp(-1i*k_u(m)*a*sin(theta));

zeta_p = fft(zeta);
psi_p = fft(psi);

A_p = k_u(m)*diff_bessel(2,p,1,k_u(m)*a)./besselh(p,k_u(m)*a); 
B_p = tau2*k_w(m)*diff_bessel(1,p,1,k_w(m)*a)./besselj(p,k_w(m)*a); 
coef_p = 1./(B_p-A_p);
Uhat = coef_p.*(B_p.*zeta_p - psi_p);
What = coef_p.*(A_p.*zeta_p - psi_p);

   U_exact(:,m) = ifft(Uhat);
   BU_exact(:,m) = ifft(Uhat.*besselh(p,k_u(m).*b)./besselh(p,k_u(m).*a));
   W_exact(:,m) = ifft(What);

U_norm(m) = norm(U_exact(:,m),2)/sqrt(N_theta);
BU_norm(m) = norm(BU_exact(:,m),2)/sqrt(N_theta);
W_norm(m) = norm(W_exact(:,m),2)/sqrt(N_theta);

end

filename = sprintf('data_%s.mat','NoPer');
save(filename,'M','lambda','U_norm','W_norm','BU_norm',...
    'lambda_crit','epsilon_u_plot','epsilon_w_plot','Eps','a','b',...
    'N_theta','N')

% lambda_crit_plot = lambda_crit*ones(M,1);
% norm_max = max([max(U_norm),max(W_norm),max(BU_norm)]);
% yy = linspace(0,norm_max,M)';


% figure(1);
% plot(lambda,U_norm,'b-o',lambda,W_norm,'g-*',...
%     lambda,BU_norm,'c-d',lambda_crit_plot,yy,'r--');
% xlabel('$\lambda$','interpreter','latex');
% ylabel('$|U|_2$ and $|W|_2$','interpreter','latex');
% title('$|U|_2$ and $|W|_2$ versus $\lambda$','interpreter','latex');
% legend('|U|_2','|W|_2','|BU|_2','lambda_c');
% 


%
% Plot of critical lambda
%

% eps_max = max(max(epsilon_u_plot),max(-real(epsilon_w_plot)));
% yy = linspace(0,eps_max,M)';
% 
% figure(2);
% plot(lambda,epsilon_u_plot,'b-o',lambda,-real(epsilon_w_plot),'g-*',...
%     lambda_crit_plot,yy,'r--');
% xlabel('$\lambda$','interpreter','latex');
% ylabel('$\epsilon$','interpreter','latex');
% title('$\epsilon_u$ and -Real($\epsilon_w$) versus $\lambda$','interpreter','latex');
% legend('epsilon_u','-Re[epsilon_w]','lambda_c');

 

