close all

L = 2*pi; 
N = 0;
N_theta = 64;
theta = (L/N_theta)*[0:N_theta-1]';
Y_p = 1i*3.4.*ones(N_theta,1); Z_p = -1i*3.4.*ones(N_theta,1);

N_lambda = 121; lambda_low = 34; lambda_high = 39;
g_bar = 1; b = 10*g_bar; Eps = 0;

lambda = linspace(lambda_low,lambda_high,N_lambda);
k_0 = (2*pi./lambda).';
n_u = zeros(N_lambda,1); n_w = zeros(N_lambda,1);
epsilon_u = zeros(N_lambda,1); epsilon_w = zeros(N_lambda,1);
sigma_hat = zeros(N_lambda,1);


%% indices
omega = k_0*(3e8)*(1e6); mu = 0.4; 
OUT = 'VACUUM'; IN = 'VACUUM';

for i = 1:N_lambda
    [n_u(i),epsilon_u(i)] = ri_perm(lambda(i),OUT);
    [n_w(i),epsilon_w(i)] = ri_perm(lambda(i),IN);
    [~,~,~,~,sigma_hat(i)] = sigma_hat_graphene_low(omega(i),300,mu,1,0.45,2.6*(1e-3));
end

k_u = n_u.*k_0; k_w = n_w.*k_0; 
sigma_u = 1./epsilon_u; sigma_w = 1./epsilon_w;
eta = sigma_hat./(1i.*k_0.*epsilon_w);

%% set up
p = [0:N_theta/2-1,-N_theta/2:-1]';
f = cos(4*theta);
f_theta = ifft(1i*p.*fft(f));
P_0 = 1; % coefficient of incident

I_u_n= zeros(N_theta,N+1,N_lambda); I_w_n= zeros(N_theta,N+1,N_lambda);
Iu_norm = zeros(N_lambda,1); Iw_norm = zeros(N_lambda,1);
Iu = zeros(N_theta,N_lambda);Iw = zeros(N_theta,N_lambda);

Q_u_n= zeros(N_theta,N+1,N_lambda); S_w_n= zeros(N_theta,N+1,N_lambda);
Qu = zeros(N_theta,N_lambda);Sw = zeros(N_theta,N_lambda);
Qu_norm = zeros(N_lambda,1); Sw_norm = zeros(N_lambda,1);


%% main 
for m=1:N_lambda
zeta_n = zeros(N_theta,N+1); psi_n = zeros(N_theta,N+1); 
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);
Sin = sin(theta); Exp = exp(-1i*k_u(m)*g_bar.*Sin);
zeta_n(:,0+1) = -P_0.*Exp; 
psi_n(:,0+1) = (1i*k_u(m))*g_bar*Sin.*Exp;

for n=1:N
    f_n = f.*f_n/n;
    if n > 1
        f_nmo = f.*f_nmo/(n-1);
    end
    zeta_n(:,n+1) = -Exp.*(-1i*k_u(m))^n.*f_n.*Sin.^n;
    psi_n(:,n+1) = (1i*k_u(m)).*Exp.*(g_bar*(-1i*k_u(m))^n.*f_n.*...
    Sin.^(n+1)+(f.*Sin-f_theta.*cos(theta)).*...
    (-1i*k_u(m))^(n-1).*f_nmo.*Sin.^(n-1));
end

Nzeta_n = current_n(f_theta,zeta_n,N);

[I_u_n(:,:,m),I_w_n(:,:,m)] = IIO_fe_TM_twolayer(Nzeta_n,psi_n,f,f_theta,p,...
    k_u(m),k_w(m),sigma_u(m),sigma_w(m),eta(m),g_bar,N_theta,N,Y_p,Z_p);
anp = field_fe_IIO_helmholtz_polar_exterior(I_u_n(:,:,m),f,f_theta,k_u(m),g_bar,p,N_theta,N,sigma_u(m),Y_p);
Q_u_n(:,:,m) = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u(m),g_bar,p,N_theta,N,sigma_u(m),Z_p);
dnp = field_fe_IIO_helmholtz_polar_interior(I_w_n(:,:,m),f,f_theta,k_w(m),g_bar,p,N_theta,N,sigma_w(m),Z_p);
S_w_n(:,:,m) = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w(m),g_bar,p,N_theta,N,sigma_w(m),Y_p);

for j=1:N_theta
    k = floor(N/2);
    Qu(j,m) = padesum(Q_u_n(j,:,m).',Eps,k);
    Sw(j,m) = padesum(S_w_n(j,:,m).',Eps,k);
    Iu(j,m) = padesum(I_u_n(j,:).',Eps,k);
    Iw(j,m) = padesum(I_w_n(j,:).',Eps,k);
end

Iu_norm(m,1) = norm(Iu(:,m),2)/sqrt(N_theta);
Iw_norm(m,1) = norm(Iw(:,m),2)/sqrt(N_theta);
Qu_norm(m,1) = norm(Qu(:,m),2)/sqrt(N_theta);
Sw_norm(m,1) = norm(Sw(:,m),2)/sqrt(N_theta);
     
end

[Qu_flat,Qu_index_flat] = max(Qu_norm(:,1));[Sw_flat,Sw_index_flat] = max(Sw_norm(:,1));
fprintf('Qu_index = %g, Sw_index = %g\n', Qu_index_flat,Sw_index_flat);

plot(lambda,Qu_norm,'b-o',lambda,Sw_norm,'r-*');
hold on
plot(lambda(Qu_index_flat),Qu_flat,'k^',lambda(Sw_index_flat),Sw_flat,'kv');
xlabel('$\lambda$','interpreter','latex');
ylabel('Qu & Sw');
title('Qu Sw versus $\lambda$','interpreter','latex');
ll = legend('Qu','Sw');
set(ll,'FontSize',16,'interpreter','latex','Location','best');

