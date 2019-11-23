function [U_norm,W_norm,BU_norm] = FE_app_TM_pade(f,N_theta,theta,...
    a,b,N,M,Eps_max,N_eps,OUT,IN,mu,lambda_low,lambda_high,name)
% this function computes and saves the data to file
%% example
% N_theta = 64;N = 16;M = 201;N_eps = 101;
% a = 0.025;b = 10*a;Eps_max = 0.01*a;
% L = 2*pi;theta = (L/N_theta)*[0:N_theta-1]';
% f = exp(cos(theta));name = 'expcos100_eps';
% OUT = 'VACUUM'; IN = 'VACUUM';
% FE_app_pade(M,f,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

lambda = linspace(lambda_low,lambda_high,M);
epsvec = linspace(0,Eps_max,N_eps);
k_0 = (2*pi./lambda).';
n_u = zeros(M,1);
n_w = zeros(M,1);
epsilon_u = zeros(M,1);
epsilon_w = zeros(M,1);
epsilon_u_plot = zeros(M,1);
epsilon_w_plot = zeros(M,1);
sigma_hat = zeros(M,1);

%% indices
omega = k_0*(3e8)*(1e6);
for i = 1:M
    [n_u(i),epsilon_u(i)] = ri_perm(lambda(i),OUT);
    [n_w(i),epsilon_w(i)] = ri_perm(lambda(i),IN);
    epsilon_u_plot(i) = epsilon_u(i);
    epsilon_w_plot(i) = epsilon_w(i);
    [sigma_hat(i),~,~,~,~] = sigma_hat_graphene_low(omega(i),300,mu,1,0.45,2.6*(1e-3));
end

% [mini,r] = min(abs(epsilon_u - real(-epsilon_w)));
% lambda_crit = lambda(r);

k_u = n_u.*k_0; 
k_w = n_w.*k_0; 
tau2 = epsilon_u./epsilon_w;
eta = sigma_hat./(1i.*k_0.*epsilon_w);

%% set up
p = [0:N_theta/2-1,-N_theta/2:-1]';
f_theta = ifft(1i*p.*fft(f));

P_0 = 1; % coefficient of incident


U_n = zeros(N_theta,N+1,M); W_n = zeros(N_theta,N+1,M); %thata and order n 
U = zeros(N_theta,M); W = zeros(N_theta,M); % theta and sum on n 
U_norm = zeros(N_eps,M); W_norm = zeros(N_eps,M);

Gn_U_n = zeros(N_theta,N+1,M); Gn_W_n = zeros(N_theta,N+1,M);
Gn_U = zeros(N_theta,M); Gn_W = zeros(N_theta,M);
Gn_U_norm = zeros(N_eps,M); Gn_W_norm = zeros(N_eps,M);

B_far = zeros(N_theta,N+1,M); BU = zeros(N_theta,M); 
BU_norm = zeros(N_eps,M);


%% main 
for m=1:M
zeta_n = zeros(N_theta,N+1); psi_n = zeros(N_theta,N+1); 
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);
Sin = sin(theta);
Exp = exp(-1i*k_u(m)*a.*Sin);
zeta_n(:,0+1) = -P_0.*Exp; 
psi_n(:,0+1) = (1i*k_u(m))*a*Sin.*Exp;

for n=1:N
    f_n = f.*f_n/n;
    if n > 1
        f_nmo = f.*f_nmo/(n-1);
    end
    zeta_n(:,n+1) = -Exp.*(-1i*k_u(m))^n.*f_n.*Sin.^n;
    psi_n(:,n+1) = (1i*k_u(m)).*Exp.*(a*(-1i*k_u(m))^n.*f_n.*...
    Sin.^(n+1)+(f.*Sin-f_theta.*cos(theta)).*...
    (-1i*k_u(m))^(n-1).*f_nmo.*Sin.^(n-1));
end

Nzeta_n = current_n(f_theta,zeta_n,N);

[U_n(:,:,m),W_n(:,:,m)] = dno_fe_TM_twolayer(Nzeta_n,psi_n,f,f_theta,p,k_u(m),k_w(m),eta(m),tau2(m),a,N_theta,N);
apn_fe = field_fe_helmholtz_polar_exterior(U_n(:,:,m),f,k_u(m),a,p,N_theta,N);
Gn_U_n(:,:,m) = dno_fe_helmholtz_polar_exterior(apn_fe,f,f_theta,k_u(m),a,p,N_theta,N);
dpn_fe = field_fe_helmholtz_polar_interior(W_n(:,:,m),f,k_w(m),a,p,N_theta,N);
Gn_W_n(:,:,m) = dno_fe_helmholtz_polar_interior(dpn_fe,f,f_theta,k_w(m),a,p,N_theta,N);   
B_far(:,:,m) = B_fe_farfield(apn_fe,k_u(m),a,b,p,N_theta,N);

   for l=1:N_eps
       Eps = epsvec(l);
       for j=1:N_theta
           k = floor(N/2);
           %BU(j,m) = padesum(B_far(j,:).',Eps,k);
           BU(j,m) = padesum(B_far(j,:,m).',Eps,k);
           U(j,m) = padesum(U_n(j,:,m).',Eps,k);
           W(j,m) = padesum(W_n(j,:,m).',Eps,k);
           Gn_U(j,m) = padesum(Gn_U_n(j,:,m).',Eps,k);
           Gn_W(j,m) = padesum(Gn_W_n(j,:,m).',Eps,k);          
       end
       
       U_norm(l,m) = norm(U(:,m),2)/sqrt(N_theta);
       BU_norm(l,m) = norm(BU(:,m),2)/sqrt(N_theta);
       W_norm(l,m) = norm(W(:,m),2)/sqrt(N_theta);
       Gn_U_norm(l,m) = norm(Gn_U(:,m),2)/sqrt(N_theta);
       Gn_W_norm(l,m) = norm(Gn_W(:,m),2)/sqrt(N_theta);
   end
end

filename = sprintf('data_%s.mat',name);
save(filename,'M','lambda','N_theta','N','N_eps','Eps_max','OUT','IN',...
    'epsilon_u_plot','epsilon_w_plot','epsvec','a','b',...
    'U_n','W_n','Gn_U_n','Gn_W_n','B_far',...
    'U_norm','W_norm','Gn_U_norm','Gn_W_norm','BU_norm')







 

