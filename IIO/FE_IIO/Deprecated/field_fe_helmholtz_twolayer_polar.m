function [anp] = field_fe_helmholtz_twolayer_polar(zeta_n,psi_n,f,f_theta,tau2,...
    p,k_u,k_w,a,N_theta,N)
% function [anp,dnp] = field_fe_helmholtz_twolayer_polar(zeta,psi,f,f_theta,tau2,...
%     p,k,a,N_theta,N)
% can simplify some matrix such as zetaDnp
%check the order of the input

%% set up
anp = zeros(N_theta,N+1); 
% dnp = zeros(N_theta,N+1); 
% fn = zeros(N_theta,N+1);
% fn(:,0+1) = ones(N_theta,1);
% for n=1:N
%   fn(:,n+1) = f.*fn(:,n-1+1)/n;
% end

UAnp = zeros(N_theta,N+1,N+1); %fe_exterior for U_0,...,U_n
UDnp = zeros(N_theta,N+1,N+1); %fe_interior for U_0,...,U_n
%W = zeros(N_theta,N+1,N+1); %fe for W_0,...,W_n
zetaDnp = zeros(N_theta,N+1,N+1); %fe for zeta_0,...,zeta_n
G_u = zeros(N_theta,N+1,N+1); %dno for G_u_jU_i
G_w = zeros(N_theta,N+1,N+1); %dno for G_w_jU_i
% zeta_u = zeros(N_theta,N+1,N+1); %dno for G_u_jzeta_i
zeta_w = zeros(N_theta,N+1,N+1); %dno for G_w_jzeta_i

%% compute order 0
Delta1 = -a*k_u*diff_besselh(p,1,a*k_u)./besselh(p,a*k_u);
Delta2 = a*k_w*diff_besselj(p,1,a*k_w)./besselj(p,a*k_w);
Delta = 1./(Delta1 + tau2*Delta2);
Qhat = fft( zeta_n(:,0+1) );
Rhat = -fft( psi_n(:,0+1) );
anp(:,0+1) = Delta.*(Rhat + tau2.*Delta2.*Qhat);
% dnp(:,0+1) = Delta.*(-Rhat - Delta1.*Qhat);
UAnp(:,:,0+1) = field_fe_helmholtz_polar_exterior(anp(:,0+1),f,k_u,p,N_theta,N,a);
UDnp(:,:,0+1) = field_fe_helmholtz_polar_interior(anp(:,0+1),f,k_w,p,N_theta,N,a);
% zetaAnp(:,:,0+1) = field_fe_helmholtz_polar_exterior(zeta_n(:,0+1),f,k_u,p,N_theta,N,a);
zetaDnp(:,:,0+1) = field_fe_helmholtz_polar_interior(zeta_n(:,0+1),f,k_w,p,N_theta,N,a);
G_u(:,1:N,0+1) = dno_fe_helmholtz_polar_exterior(UAnp(:,:,0+1),f,f_theta,p,a,k_u,N_theta,N-1);
G_w(:,1:N,0+1) = dno_fe_helmholtz_polar_interior(UDnp(:,:,0+1),f,f_theta,p,a,k_w,N_theta,N-1);
% zeta_u(:,0+1,0+1) = dno_fe_helmholtz_polar_exterior(zetaAnp(:,:,0+1),f,f_theta,p,a,k_u,N_theta,0+1);
zeta_w(:,:,0+1) = dno_fe_helmholtz_polar_interior(zetaDnp(:,:,0+1),f,f_theta,p,a,k_w,N_theta,N);


%% order N
for n=1:N
  Rhat = -fft( psi_n(:,n+1) );
  Qhat = 0;
  zetaDnp(:,:,n+1)=field_fe_helmholtz_polar_interior(zeta_n(:,n+1),f,k_w,p,N_theta,N,a);
  zeta_w(:,1:N+1-n,n+1) = dno_fe_helmholtz_polar_interior(zetaDnp(:,:,n+1),f,f_theta,p,a,k_w,N_theta,N-n);
  for m=0:n
    Rhat = Rhat + tau2*fft(zeta_w(:,n+1-m,m+1));
  end
  for m=0:n-1
    Qhat = fft(G_u(:,n+1-m,m+1)) + tau2*fft(G_w(:,n+1-m,m+1));
  end
  anp(:,n+1) = Delta.*(Rhat + Qhat);
  UAnp(:,:,n+1) = field_fe_helmholtz_polar_exterior(anp(:,n+1),f,k_u,p,N_theta,N,a);
  UDnp(:,:,n+1) = field_fe_helmholtz_polar_interior(anp(:,n+1),f,k_w,p,N_theta,N,a);
  G_u(:,1:N+1-n,n+1) = dno_fe_helmholtz_polar_exterior(UAnp(:,:,n+1),f,f_theta,p,a,k_u,N_theta,N-n);
  G_w(:,1:N+1-n,n+1) = dno_fe_helmholtz_polar_interior(UDnp(:,:,n+1),f,f_theta,p,a,k_w,N_theta,N-n);
end