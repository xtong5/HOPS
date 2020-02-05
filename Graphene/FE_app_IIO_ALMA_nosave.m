function [Qu_norm,Sw_norm,BU_norm] = FE_app_IIO_ALMA_nosave(f,N_theta,theta,...
    a,b,Y_p,Z_p,N,Eps_max,N_eps,OUT,IN,mu,lambda_low,lambda_high,N_lambda)
% this function computes and saves the data to file
%% example
% N_theta = 64;N = 16;M = 201;N_eps = 101;
% a = 0.025;b = 10*a;Eps_max = 0.01*a;
% L = 2*pi;theta = (L/N_theta)*[0:N_theta-1]';
% f = exp(cos(theta));name = 'expcos100_eps';
% OUT = 'VACUUM'; IN = 'VACUUM';
% FE_app_pade(M,f,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

lambda = linspace(lambda_low,lambda_high,N_lambda);
epsvec = linspace(0,Eps_max,N_eps);
k_0 = (2*pi./lambda).';
n_u = zeros(N_lambda,1); n_w = zeros(N_lambda,1);
epsilon_u = zeros(N_lambda,1); epsilon_w = zeros(N_lambda,1);
% epsilon_u_plot = zeros(N_lambda,1); epsilon_w_plot = zeros(N_lambda,1);
sigma_hat = zeros(N_lambda,1);

%% indices
omega = k_0*(3e8)*(1e6);
for i = 1:N_lambda
    [n_u(i),epsilon_u(i)] = ri_perm(lambda(i),OUT);
    [n_w(i),epsilon_w(i)] = ri_perm(lambda(i),IN);
%     epsilon_u_plot(i) = epsilon_u(i); epsilon_w_plot(i) = epsilon_w(i);
    [sigma_hat(i),~,~,~,~] = sigma_hat_graphene_low(omega(i),300,mu,1,0.45,2.6*(1e-3));
end

% [mini,r] = min(abs(epsilon_u - real(-epsilon_w)));
% lambda_crit = lambda(r);

k_u = n_u.*k_0; k_w = n_w.*k_0; 
sigma_u = 1./epsilon_u; sigma_w = 1./epsilon_w;
eta = sigma_hat./(1i.*k_0.*epsilon_w);

%% set up
p = [0:N_theta/2-1,-N_theta/2:-1]';
f_theta = ifft(1i*p.*fft(f));

P_0 = 1; % coefficient of incident

I_u_n= zeros(N_theta,N+1,N_lambda); I_w_n= zeros(N_theta,N+1,N_lambda);
Iu_norm = zeros(N_eps,N_lambda); Iw_norm = zeros(N_eps,N_lambda);
Iu = zeros(N_theta,N_lambda);Iw = zeros(N_theta,N_lambda);

Q_u_n= zeros(N_theta,N+1,N_lambda); S_w_n= zeros(N_theta,N+1,N_lambda);
Qu = zeros(N_theta,N_lambda);Sw = zeros(N_theta,N_lambda);
Qu_norm = zeros(N_eps,N_lambda); Sw_norm = zeros(N_eps,N_lambda);

B_far = zeros(N_theta,N+1,N_lambda); BU = zeros(N_theta,N_lambda); 
BU_norm = zeros(N_eps,N_lambda);

%% main 
for m=1:N_lambda
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

[I_u_n(:,:,m),I_w_n(:,:,m)] = IIO_fe_TM_twolayer(Nzeta_n,psi_n,f,f_theta,p,...
    k_u(m),k_w(m),sigma_u(m),sigma_w(m),eta(m),a,N_theta,N,Y_p,Z_p);
anp = field_fe_IIO_helmholtz_polar_exterior(I_u_n(:,:,m),f,f_theta,k_u(m),a,p,N_theta,N,sigma_u(m),Y_p);
Q_u_n(:,:,m) = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u(m),a,p,N_theta,N,sigma_u(m),Z_p);
dnp = field_fe_IIO_helmholtz_polar_interior(I_w_n(:,:,m),f,f_theta,k_w(m),a,p,N_theta,N,sigma_w(m),Z_p);
S_w_n(:,:,m) = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w(m),a,p,N_theta,N,sigma_w(m),Y_p);
B_far(:,:,m) = B_fe_farfield(anp,k_u(m),a,b,p,N_theta,N);

   for l=1:N_eps
       Eps = epsvec(l);
       for j=1:N_theta
           k = floor(N/2);
           BU(j,m) = padesum(B_far(j,:,m).',Eps,k);
           Qu(j,m) = padesum(Q_u_n(j,:,m).',Eps,k);
           Sw(j,m) = padesum(S_w_n(j,:,m).',Eps,k);
           Iu(j,m) = padesum(I_u_n(j,:).',Eps,k);
           Iw(j,m) = padesum(I_w_n(j,:).',Eps,k); 
        end
       Iu_norm(l,m) = norm(Iu(:,m),2)/sqrt(N_theta);
       BU_norm(l,m) = norm(BU(:,m),2)/sqrt(N_theta);
       Iw_norm(l,m) = norm(Iw(:,m),2)/sqrt(N_theta);
       Qu_norm(l,m) = norm(Qu(:,m),2)/sqrt(N_theta);
       Sw_norm(l,m) = norm(Sw(:,m),2)/sqrt(N_theta);
   end
end

