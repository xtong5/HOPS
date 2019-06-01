function [lambda_crit,U_norm,W_norm,BU_norm] = ...
    FE_applicationTaylor(M,f,N_theta,a,b,N,Eps,theta,name)
% this function computes and saves the data to file
% data given by different lambda (taylorsum)


lambda = linspace(0.3,0.8,M);
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

[mini,r] = min(abs(epsilon_u - real(-epsilon_w)));
lambda_crit = lambda(r);

k_u = n_u*k_zero; 
k_w = n_w.*k_zero; 
Mode = 2; 
% Mode = 1;

fprintf('test\n');
fprintf('-------------\n');
fprintf('Eps = %g  a = %g  b = %g\n',Eps,a,b);
fprintf('N_theta = %d N = %d M = %d\n',N_theta,N,M);
fprintf('\n');

p = [0:N_theta/2-1,-N_theta/2:-1]';
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

% zeta = -P_0*exp(-1i*k_u(m).*A.*sin(theta));
% psi = (1i*k_u(m)).*A.*zeta;
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
     BU(j,m) = taylorsum(B_far(j,:).',Eps,N);
     U(j,m) = taylorsum(U_n(j,:,m).',Eps,N);
     W(j,m) = taylorsum(W_n(j,:,m).',Eps,N);
   end
   
   U_norm(m) = norm(U(:,m),2)/sqrt(N_theta);
   BU_norm(m) = norm(BU(:,m),2)/sqrt(N_theta);
   W_norm(m) = norm(W(:,m),2)/sqrt(N_theta);
   
end

filename = sprintf('data_%s.mat',name);
save(filename,'M','lambda','U_norm','W_norm','BU_norm',...
    'lambda_crit','epsilon_u_plot','epsilon_w_plot','Eps','a','b',...
    'N_theta','N')




 

