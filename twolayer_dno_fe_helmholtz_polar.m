function [U_n] = twolayer_dno_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,tau2,...
    p,k_u,k_w,a,N_theta,N)
% function [anp,dnp] = field_fe_helmholtz_twolayer_polar(zeta,psi,f,f_theta,tau2,...
%     p,k,a,N_theta,N)
% can simplify some matrix such as zetaDnp
%check the order of the input

%% set up
U_n = zeros(N_theta,N+1);
G_u_U = zeros(N_theta,N+1,N+1); %dno for G_u_jU_i
G_w_U = zeros(N_theta,N+1,N+1); %dno for G_w_jU_i
G_w_zeta = zeros(N_theta,N+1,N+1); %dno for G_w_jzeta_i
Delta1 = -a*k_u*diff_besselh(p,1,a*k_u)./besselh(p,a*k_u);
Delta2 = a*k_w*diff_besselj(p,1,a*k_w)./besselj(p,a*k_w);
Delta = 1./(Delta1 + tau2*Delta2);
Uanp = zeros(N_theta,N+1,N+1);
Udnp = zeros(N_theta,N+1,N+1);
zetadnp = zeros(N_theta,N+1,N+1);


%% Compute and store G_w[zeta]
for n=0:N
  s = n;
  xi = zeros(N_theta,N-s+1);
  xi(:,0+1) = zeta_n(:,s+1);
%   zetadnp(:,1:N-s+1,s+1) = field_fe_helmholtz_polar_interior(xi,f,k_w,a,p,N_theta,N-s);
%   G_w_zeta(:,1:N-s+1,s+1) = dno_fe_helmholtz_polar_interior(zetadnp(:,:,s+1),f,f_theta,k_w,a,p,N_theta,N-s);
  dnp = field_fe_helmholtz_polar_interior(xi,f,k_w,a,p,N_theta,N-s);
  G_w_zeta(:,1:N-s+1,s+1) = dno_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N-s);
end


%% compute order n=0
Qhat = fft( zeta_n(:,0+1) );
Rhat = -fft( psi_n(:,0+1) );
U_n(:,0+1) = ifft(Delta.*(Rhat + tau2.*Delta2.*Qhat));
% compute and store G_u[U_0] and G_w[U_0]
xi = zeros(N_theta,N+1);
xi(:,0+1) = U_n(:,0+1);
% Uanp(:,1:N,0+1) = field_fe_helmholtz_polar_exterior(xi,f,k_u,a,p,N_theta,N-1);
% G_u_U(:,1:N,0+1) = dno_fe_helmholtz_polar_exterior(Uanp(:,:,0+1),f,f_theta,k_u,a,p,N_theta,N-1);
% Udnp(:,1:N,0+1) = field_fe_helmholtz_polar_interior(xi,f,k_w,a,p,N_theta,N-1);
% G_w_U(:,1:N,0+1) = dno_fe_helmholtz_polar_interior(Udnp(:,:,0+1),f,f_theta,k_w,a,p,N_theta,N-1);
anp = field_fe_helmholtz_polar_exterior(xi,f,k_u,a,p,N_theta,N-1);
G_u_U(:,1:N,0+1) = dno_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N-1);
dnp = field_fe_helmholtz_polar_interior(xi,f,k_w,a,p,N_theta,N-1);
G_w_U(:,1:N,0+1) = dno_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N-1);


%% order n>0
for n=1:N
  s=n;
  Rhat = -fft( psi_n(:,n+1) );
  Qhat = 0;
  for m=0:n
    Rhat = Rhat + tau2*fft(G_w_zeta(:,n+1-m,m+1));
  end
  for m=0:n-1
    Qhat = fft(G_u_U(:,n+1-m,m+1)) + tau2*fft(G_w_U(:,n+1-m,m+1));
  end
  U_n(:,n+1) = ifft(Delta.*(Rhat + Qhat));
  % compute and store G_u[U_s] and G_w[U_s]
  xi = zeros(N_theta,N-s+1);
  xi(:,0+1) = U_n(:,s+1);
%   Uanp(:,1:N-s+1,s+1) = field_fe_helmholtz_polar_exterior(xi,f,k_u,a,p,N_theta,N-s);
%   Udnp(:,1:N-s+1,s+1) = field_fe_helmholtz_polar_interior(xi,f,k_w,a,p,N_theta,N-s);
%   G_u_U(:,1:N-s+1,s+1) = dno_fe_helmholtz_polar_exterior(Uanp(:,:,s+1),f,f_theta,k_u,a,p,N_theta,N-s);
%   G_w_U(:,1:N-s+1,s+1) = dno_fe_helmholtz_polar_interior(Udnp(:,:,s+1),f,f_theta,k_w,a,p,N_theta,N-s);
  anp = field_fe_helmholtz_polar_exterior(xi,f,k_u,a,p,N_theta,N-s);
  G_u_U(:,1:N-s+1,s+1) = dno_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N-s);
  dnp = field_fe_helmholtz_polar_interior(xi,f,k_w,a,p,N_theta,N-s);
  G_w_U(:,1:N-s+1,s+1) = dno_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N-s); 
end