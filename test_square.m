clear all;
SavePlots = 0;

SumType = 1;

warning('off')

L = 2*pi;
N_theta = 64;
a = 0.025;
b = 10*a;
N = 16;
Eps = 1; %0.1*a
M=1;

Mode = 2; %check 
% Mode = 1;

%N_theta_vec=[16, 64, 128, 256];
N_theta_vec=[28*8, 29*8, 30*8];
for l=1:length(N_theta_vec)
    N_theta = N_theta_vec(l);
theta = (L/N_theta)*[0:N_theta-1]';

ff = zeros(N_theta/4,1);
for i=0:N_theta/8
    ff(i+1) = a/cos(theta(i+1));
end
for i=N_theta/8+1:N_theta/8*2-1
    ff(i+1) = a/sin(theta(i+1));
end
f = [ff; ff; ff; ff];

lambda = 0.45;
k_zero = 2*pi./lambda;
n_u = 2.5;
m_w = zeros(M,1);
epsilon_w = zeros(M,1);
epsilon_u_plot = zeros(M,1);
epsilon_w_plot = zeros(M,1);
epsilon_u = n_u^2;

[n_w,epsilon_w] = ri_perm(lambda,'SILVER');
epsilon_u_plot = epsilon_u;
epsilon_w_plot = epsilon_w;

lambda_crit = lambda;

k_u = n_u*k_zero; 
k_w = n_w.*k_zero; 

fprintf('test\n');
fprintf('-------------\n');
fprintf('Eps = %g  a = %g  b = %g\n',Eps,a,b);
fprintf('lambda = %g  k_u = %g  k_w = %g\n',lambda,k_u,k_w);
fprintf('N_theta = %d N = %d M = %d\n',N_theta,N,M);
fprintf('\n');

p = [0:N_theta/2-1,-N_theta/2:-1]';
f_theta = ifft(1i*p.*fft(f));


P_0 = 1;

BU = zeros(N_theta,1);
U = zeros(N_theta,1);
W = zeros(N_theta,1);
B_far = zeros(N_theta,N+1);
% BU_norm = zeros(M,1);
% U_norm = zeros(M,1);
% W_norm = zeros(M,1);
U_n = zeros(N_theta,N+1);
W_n = zeros(N_theta,N+1); 


% zeta = -P_0*exp(-1i*k_u.*A.*sin(theta));
% psi = (1i*k_u).*A.*zeta;
zeta_n = zeros(N_theta,N+1);
psi_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);
Sin = sin(theta);
Exp = exp(-1i*k_u*a.*Sin);
zeta_n(:,0+1) = -Exp; 
psi_n(:,0+1) = (1i*k_u)*a*Sin.*Exp;
for n=1:N
f_n = f.*f_n/n;
if n > 1
   f_nmo = f.*f_nmo/(n-1);
end

zeta_n(:,n+1) = -Exp.*(-1i*k_u)^n.*f_n.*Sin.^n;
psi_n(:,n+1) = (1i*k_u).*Exp.*(a*(-1i*k_u)^n.*f_n.*Sin.^(n+1)+...
    (f.*Sin-f_theta.*cos(theta)).*(-1i*k_u)^(n-1).*f_nmo.*Sin.^(n-1));
end

if(Mode==1)
  tau2 = 1;
else
  tau2 = k_u^2/k_w^2;
end


   U_n(:,:) = twolayer_dno_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,tau2,...
    p,k_u,k_w,a,N_theta,N);
   apn_fe = field_fe_helmholtz_polar_exterior(U_n(:,:),f,k_u,a,p,N_theta,N);
   W_n(:,:) = U_n(:,:) - zeta_n;
   
   for n=1:N+1
     B_far(:,n)=ifft(apn_fe(:,n).*besselh(p,k_u.*b)./besselh(p,k_u.*a));
   end 

  for j=1:N_theta
     k = floor(N/2);
     if(SumType==1)
         BU(j,1) = taylorsum(B_far(j,:).',Eps,N);
         U(j,1) = taylorsum(U_n(j,:).',Eps,N);
         W(j,1) = taylorsum(W_n(j,:).',Eps,N);
     else
         BU(j,1) = padesum(B_far(j,:).',Eps,k);
         U(j,1) = padesum(U_n(j,:).',Eps,k);
         W(j,1) = padesum(W_n(j,:).',Eps,k);
     end
   end
   
   
   U_norm(l) = norm(U(:),2)/sqrt(N_theta);
   BU_norm(l) = norm(BU(:),2)/sqrt(N_theta);
   W_norm(l) = norm(W(:),2)/sqrt(N_theta);

end

for j=1:length(N_theta_vec)
fprintf('U_norm = %20.16g, W_norm = %20.16g, BU_norm = %20.16g \n',...
    U_norm(j),W_norm(j),BU_norm(j));
end