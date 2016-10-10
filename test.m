clear all
close all

L = 2*pi;
%lambda = [0.2:0.005:1];
lambda = linspace(0.2,1.0,11);
M = size(lambda,2);
% lambda = 0.2;
k_zero = 2*pi./lambda;
n_u = 1;
% n_w = 1.3;
for i = 1:M
   n_w(i) = ri_perm(lambda(i),'SILVER');
end
k_u = n_u*k_zero; 
k_w = n_w.*k_zero; 
RunNumber = 0;
%Mode = 2; %check 
Mode = 1;

if(RunNumber==0)
  Eps = 0;
  N_theta = 8;
  a = 1;
  b = 5;
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

f = exp(cos(theta));
f_theta = -sin(theta).*f;
A = a+Eps*f;

BU = zeros(N_theta,M);
U = zeros(N_theta,M);
for m =1:M
% u_i = exp(-1i*k_u(m).*A);
% zeta = -u_i;
% psi = A.*(1i*k_u(m)).*exp(-1i*k_u(m).*A);
zeta_n = zeros(N_theta,N+1);
psi_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);
Sin = sin(theta);
zeta_n(:,0+1) = exp(-1i*k_u(m)*a).*f_n; 
psi_n(:,0+1) = (1i*k_u(m))*a*sin(theta).*exp(-1i*k_u(m)*a.*Sin);
for n=1:N
f_n = f.*f_n/n;
if n > 1
   f_nmo = f.*f_nmo/(n-1);
end
Exp = exp(-1i*k_u(m)*a.*Sin);
zeta_n(:,n+1) = -Exp.*(-1i*k_u(m))^n.*f_n.*Sin.^n;
psi_n(:,n+1) = (1i*k_u(m)).*Exp.*(a*(-1i*k_u(m))^n.*f_n.*Sin.^(n+1)+...
    (f.*Sin-f_theta.*cos(theta)).*(-1i*k_u(m))^(n-1).*f_nmo.*Sin.^(n-1));
end

if(Mode==1)
  tau2 = 1;
else
  tau2 = k_u(m)^2/k_w(m)^2;
end

   U_n = zeros(N_theta,N+1,M);
   W_n = zeros(N_theta,N+1,M); 
   U_n(:,:,m) = twolayer_dno_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,tau2,...
    p,k_u(m),k_w(m),a,N_theta,N);
   apn_fe = field_fe_helmholtz_polar_exterior(U_n(:,:,m),f,k_u(m),a,p,N_theta,N);
%    W_n(:,:,m) = U_n(:,:,m) - zeta_n;
   for n=1:N+1
     B_far(:,n)=ifft(apn_fe(:,n).*besselh(p,k_u(m).*b)./besselh(p,k_u(m).*a));
   end 

   for j=1:N_theta
     BU(j,m) = taylorsum(B_far(j,:),Eps,N);
     U(j,m) = taylorsum(U_n(j,:,m),Eps,N);
   end
   BU_norm(m) = norm(BU,inf);
   U_norm(m) = norm(U,inf);
end
plot(lambda,U_norm,'b-*',lambda,BU_norm,'r-o')



 

