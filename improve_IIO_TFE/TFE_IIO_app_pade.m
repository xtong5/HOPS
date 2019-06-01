function [lambda_crit,U_norm,W_norm,BU_norm] = ...
    TFE_IIO_app_pade(M,f,N_theta,theta,a,b,c,N,N_r,Eps_max,N_eps,OUT,IN,eta,name)
% this function computes and saves the data to file
%% example
% N_theta = 64;N = 16;N_r = 16;M = 201;N_eps = 101;
% a = 0.025;b = 10*a;Eps_max = 0.01*a;
% L = 2*pi;theta = (L/N_theta)*[0:N_theta-1]';
% f = exp(cos(theta));name = 'expcos100_eps';
% OUT = 'VACUUM'; IN = 'SILVER';
% FE_app_pade(M,f,N_theta,theta,a,b,N,Eps_max,N_eps,OUT,IN,name);

lambda = linspace(0.3,0.8,M);
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
f_theta = ifft(1i*p.*fft(f));

P_0 = 1;

U_n = zeros(N_theta,N+1,M); W_n = zeros(N_theta,N+1,M); %thata and order n 
U = zeros(N_theta,M); W = zeros(N_theta,M); % theta and sum on n 
U_norm = zeros(N_eps,M); W_norm = zeros(N_eps,M);

I_u_n= zeros(N_theta,N+1,M); I_w_n= zeros(N_theta,N+1,M);
Iu_norm = zeros(N_eps,M); Iw_norm = zeros(N_eps,M);
Iu = zeros(N_theta,M);Iw = zeros(N_theta,M);

Q_u_n= zeros(N_theta,N+1,M); S_w_n= zeros(N_theta,N+1,M);
Qu = zeros(N_theta,M);Sw = zeros(N_theta,M);
Qu_norm = zeros(N_eps,M); Sw_norm = zeros(N_eps,M);


B_far = zeros(N_theta,N+1,M); BU = zeros(N_theta,M); 
BU_norm = zeros(N_eps,M);

Gn_U_n = zeros(N_theta,N+1,M); Gn_W_n = zeros(N_theta,N+1,M);
Gn_U = zeros(N_theta,M); Gn_W = zeros(N_theta,M);
Gn_U_norm = zeros(N_eps,M); Gn_W_norm = zeros(N_eps,M);

for m=1:M
    if(Mode==1)
        sigma_u = 1;
        sigma_w = 1;
    else
        sigma_u = 1/epsilon_u(m);
        sigma_w = 1/epsilon_w(m);
    end

zeta_n = zeros(N_theta,N+1);
psi_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);
Sin = sin(theta);
Exp = exp(-1i*k_u(m)*a.*Sin);
zeta_n(:,0+1) = -P_0*Exp; 
psi_n(:,0+1) = P_0*(1i*k_u(m))*a*Sin.*Exp;

for n=1:N
    f_n = f.*f_n/n;
    if n > 1
        f_nmo = f.*f_nmo/(n-1);
    end
    zeta_n(:,n+1) = -P_0.*Exp.*(-1i*k_u(m))^n.*f_n.*Sin.^n;
    psi_n(:,n+1) = (P_0*1i*k_u(m)).*Exp.*(a*(-1i*k_u(m))^n.*f_n.*...
    Sin.^(n+1)+(f.*Sin-f_theta.*cos(theta)).*...
    (-1i*k_u(m))^(n-1).*f_nmo.*Sin.^(n-1));
end

% Z_p = sigma_u * k_u(m) * diff_besselh(p,1,k_u(m)*a)./besselh(p,k_u(m)*a);
% Y_p = sigma_w * k_w(m) * diff_besselj(p,1,k_w(m)*a)./besselj(p,k_w(m)*a);


Z_p = -1i*eta.*ones(N_theta,1);
Y_p = 1i*eta.*ones(N_theta,1);

[I_u_n(:,:,m),I_w_n(:,:,m)] = twolayer_IIO_TFE(zeta_n,psi_n,f,f_theta,p,k_u(m),k_w(m),sigma_u,sigma_w,a,b,c,N_theta,N,N_r,eta);
[Un(:,:,m),Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u(m),a,b,p,N_theta,N,N_r,sigma_u,eta);
Q_u_n(:,:,m) = IIO_tfe_helmholtz_polar_exterior(Un(:,:,m),Dr_Un,Dp_Un,f,f_theta,a,b,N,sigma_u,eta);
[Wn(:,:,m),Dr_Wn,Dp_Wn] = field_tfe_IIO_helmholtz_polar_interior(I_w_n,f,f_theta,k_w(m),a,c,p,N_theta,N,N_r,sigma_w,eta);
S_w_n(:,:,m) = IIO_tfe_helmholtz_polar_interior(Wn(:,:,m),Dr_Wn,Dp_Wn,f,f_theta,a,c,N,sigma_w,eta);

[I_u_n(:,:,m),I_w_n(:,:,m)] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u(m),k_w(m),sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);

apn = field_IIO_fe_helmholtz_polar_exterior(I_u_n(:,:,m),f,f_theta,k_u(m),a,p,N_theta,N,sigma_u,Y_p);
U_n(:,:,m) = IIO_fe_Un(apn,f,k_u(m),a,p,N_theta,N);
Gn_U_n(:,:,m) = dno_fe_helmholtz_polar_exterior(apn,f,f_theta,k_u(m),a,p,N_theta,N);
Q_u_n(:,:,m) = IIO_fe_helmholtz_polar_exterior(apn,f,f_theta,k_u(m),a,p,N_theta,N,sigma_u,Z_p);

W_n(:,:,m) = U_n(:,:,m) - zeta_n;
dpn = field_IIO_fe_helmholtz_polar_interior(I_w_n(:,:,m),f,f_theta,k_w(m),a,p,N_theta,N,sigma_w,Z_p);
Gn_W_n(:,:,m) = dno_fe_helmholtz_polar_interior(dpn,f,f_theta,k_w(m),a,p,N_theta,N);
S_w_n(:,:,m) = IIO_fe_helmholtz_polar_interior(dpn,f,f_theta,k_w(m),a,p,N_theta,N,sigma_w,Y_p);

B_far(:,:,m) = IIO_fe_farfield(apn,k_u(m),a,b,p,N_theta,N);

  
for l=1:N_eps
       Eps = epsvec(l);
       for j=1:N_theta
           k = floor(N/2);
           BU(j,m) = padesum(B_far(j,:,m).',Eps,k);
           U(j,m) = padesum(U_n(j,:,m).',Eps,k);
           W(j,m) = padesum(W_n(j,:,m).',Eps,k);
           Qu(j,m) = padesum(Q_u_n(j,:,m).',Eps,k);
           Sw(j,m) = padesum(S_w_n(j,:,m).',Eps,k);
           Iu(j,m) = padesum(I_u_n(j,:).',Eps,k);
           Iw(j,m) = padesum(I_w_n(j,:).',Eps,k); 
           Gn_U(j,m) = padesum(Gn_U_n(j,:,m).',Eps,k);
           Gn_W(j,m) = padesum(Gn_W_n(j,:,m).',Eps,k);   
       end
       U_norm(l,m) = norm(U(:,m),2)/sqrt(N_theta);
       Iu_norm(l,m) = norm(Iu(:,m),2)/sqrt(N_theta);
       BU_norm(l,m) = norm(BU(:,m),2)/sqrt(N_theta);
       W_norm(l,m) = norm(W(:,m),2)/sqrt(N_theta);
       Iw_norm(l,m) = norm(Iw(:,m),2)/sqrt(N_theta);
       Qu_norm(l,m) = norm(Qu(:,m),2)/sqrt(N_theta);
       Sw_norm(l,m) = norm(Sw(:,m),2)/sqrt(N_theta);
       Gn_U_norm(l,m) = norm(Gn_U(:,m),2)/sqrt(N_theta);
       Gn_W_norm(l,m) = norm(Gn_W(:,m),2)/sqrt(N_theta);
   end
end

filename = sprintf('data_%s.mat',name);
% save(filename,'M','lambda','U_norm','W_norm','Qu_norm','Sw_norm',...
%     'lambda_crit','epsilon_u_plot','epsilon_w_plot','epsvec','a','b',...
%     'N_theta','N','N_eps','Eps_max','OUT','IN')
save(filename,'M','lambda','N_theta','N','N_eps','Eps_max','OUT','IN',...
    'lambda_crit','epsilon_u_plot','epsilon_w_plot','epsvec','a','b','Y_p','Z_p',...
    'U_n','W_n','I_u_n','I_w_n','Q_u_n','S_w_n','Gn_U_n','Gn_W_n','B_far',...
    'U_norm','W_norm','Iu_norm','Iw_norm','Qu_norm','Sw_norm','Gn_U_norm','Gn_W_norm','BU_norm')








 

