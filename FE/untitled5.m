% compute and save the date to file 
L = 2*pi;
a = 0.025;
b = 10*a;
M = 1;
N_eps = 1;
N_theta = 64;
theta = (L/N_theta)*[0:N_theta-1]';

warning('off')

f8 = cos(8*theta);


%% N_theta96 N24
N_theta = 64;N = 16;
OUT = 'WATER';
IN = 'SILVER';
Eps_max = 0.1*a;
name = 'cos5';
m=1;

lambda = 0.3875;
epsvec = linspace(0,Eps_max,N_eps);
k_zero = (2*pi./lambda).';
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
Mode = 2; 
% Mode = 1;


p = [0:N_theta/2-1,-N_theta/2:-1]';
f8_theta = ifft(1i*p.*fft(f8));

P_0 = 1;

% BU = zeros(N_theta,M);
U = zeros(N_theta,M);
W = zeros(N_theta,M);
Gn_U = zeros(N_theta,M);
Gn_W = zeros(N_theta,M);
% B_far = zeros(N_theta,N+1);
% BU_norm = zeros(N_eps,M);
U_n = zeros(N_theta,N+1,M);
W_n = zeros(N_theta,N+1,M); 
Gn_U_n = zeros(N_theta,N+1,M);
Gn_W_n = zeros(N_theta,N+1,M);
U_norm = zeros(N_eps,M);
W_norm = zeros(N_eps,M);
Gn_U_norm = zeros(N_eps,M);
Gn_W_norm = zeros(N_eps,M);

zeta_n = zeros(N_theta,N+1);
psi_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);
Sin = sin(theta);
Exp = exp(-1i*k_u*a.*Sin);
zeta_n(:,0+1) = -Exp; 
psi_n(:,0+1) = (1i*k_u)*a*Sin.*Exp;

for n=1:N
    f_n = f8.*f_n/n;
    if n > 1
        f_nmo = f8.*f_nmo/(n-1);
    end
    zeta_n(:,n+1) = -Exp.*(-1i*k_u(m))^n.*f_n.*Sin.^n;
    psi_n(:,n+1) = (1i*k_u(m)).*Exp.*(a*(-1i*k_u(m))^n.*f_n.*...
    Sin.^(n+1)+(f8.*Sin-f8_theta.*cos(theta)).*...
    (-1i*k_u(m))^(n-1).*f_nmo.*Sin.^(n-1));
 end
if(Mode==1)
    tau2 = 1;
else
    tau2 = k_u(m)^2/k_w(m)^2;
end

  U_n(:,:,m) = twolayer_dno_fe_helmholtz_polar(zeta_n,psi_n,f8,f8_theta,tau2,...
    p,k_u(m),k_w(m),a,N_theta,N);
  apn_fe = field_fe_helmholtz_polar_exterior(U_n(:,:,m),f8,k_u(m),a,p,N_theta,N);
  Gn_U_n(:,:,m) = dno_fe_helmholtz_polar_exterior(apn_fe,f8,f8_theta,k_u(m),a,p,N_theta,N);
  W_n(:,:,m) = U_n(:,:,m) - zeta_n;
  dpn_fe = field_fe_helmholtz_polar_interior(W_n(:,:,m),f8,k_w(m),a,p,N_theta,N);
  Gn_W_n(:,:,m) = dno_fe_helmholtz_polar_interior(dpn_fe,f8,f8_theta,k_w(m),a,p,N_theta,N);
   
%    for n=1:N+1
%      B_far(:,n)=ifft(apn_fe(:,n).*besselh(p,k_u(m).*b)./besselh(p,k_u(m).*a));
%    end 

   for l=1:N_eps
       Eps = epsvec(l);
       for j=1:N_theta
           k = floor(N/2);
%            BU(j,m) = padesum(B_far(j,:).',Eps,k);
           U(j,m) = padesum(U_n(j,:,m).',Eps,k);
           W(j,m) = padesum(W_n(j,:,m).',Eps,k);
           Gn_U(j,m) = padesum(Gn_U_n(j,:,m).',Eps,k);
           Gn_W(j,m) = padesum(Gn_W_n(j,:,m).',Eps,k);          
       end
       
       U_norm(l,m) = norm(U(:,m),2)/sqrt(N_theta);
%        BU_norm(l,m) = norm(BU(:,m),2)/sqrt(N_theta);
       W_norm(l,m) = norm(W(:,m),2)/sqrt(N_theta);
       Gn_U_norm(l,m) = norm(Gn_U(:,m),2)/sqrt(N_theta);
       Gn_W_norm(l,m) = norm(Gn_W(:,m),2)/sqrt(N_theta);
   end
   

