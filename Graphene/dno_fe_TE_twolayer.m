function [U_n] = dno_fe_TE_twolayer(zeta_n,psi_n,f,f_theta,rho,...
    p,k_u,k_w,a,N_theta,N)

%% set up
U_n = zeros(N_theta,N+1);
Gn_u_U = zeros(N_theta,N+1,N+1); %dno for G_u_jU_i
Gn_w_U = zeros(N_theta,N+1,N+1); %dno for G_w_jU_i
Gn_w_zeta = zeros(N_theta,N+1,N+1); %dno for G_w_jzeta_i
Nn_U = zeros(N_theta,N+1); %normal operator for N_jU_i
Nn_zeta = zeros(N_theta,N+1); %normal operator for N_jzeta_i

Delta1 = -a*k_u*diff_besselh(p,1,a*k_u)./besselh(p,a*k_u);
Delta2 = a*k_w*diff_besselj(p,1,a*k_w)./besselj(p,a*k_w);
Delta = 1./(Delta1 + Delta2 - rho);

%% Compute and store Gn_w[zeta] and |N|n[zeta]
for n=0:N
  s = n;
  xi = zeros(N_theta,N-s+1);
  xi(:,0+1) = zeta_n(:,s+1);
  dnp = field_fe_helmholtz_polar_interior(xi,f,k_w,a,p,N_theta,N-s);
  Gn = dno_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N-s);
  Nn = rho.*current_n(f_theta,xi,N-s);
  for r=0:N-s
    Gn_w_zeta(:,r+1,s+1) = Gn(:,r+1);
    Nn_zeta(:,r+1,s+1) = Nn(:,r+1);
  end
end


%% compute order n=0
n = 0;
s = n;
Qhat = fft( zeta_n(:,0+1) );
Rhat = fft( psi_n(:,0+1) );
U_n(:,0+1) = ifft((Delta2.*Qhat - rho.*Qhat - Rhat).*Delta);

% compute and store Gn_u[U_0] and |N|n[U_0]
xi = zeros(N_theta,N-s+1);
xi(:,0+1) = U_n(:,s+1);
anp = field_fe_helmholtz_polar_exterior(xi,f,k_u,a,p,N_theta,N-s);
Gn = dno_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N-s);
Nn = rho.*current_n(f_theta,xi,N-s);
for r=0:N-s
  Gn_u_U(:,r+1,s+1) = Gn(:,r+1);
  Nn_U(:,r+1,s+1) = Nn(:,r+1);
end

% compute and store Gn_w[U_0]
xi = zeros(N_theta,N-s+1);
xi(:,0+1) = U_n(:,s+1);
dnp = field_fe_helmholtz_polar_interior(xi,f,k_w,a,p,N_theta,N-s);
Gn = dno_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N-s);
for r=0:N-s
  Gn_w_U(:,r+1,s+1) = Gn(:,r+1);
end

%% order n>0
for n=1:N
  s=n;
  Rhat = 0;
  Qhat = 0;
  for m=0:n
    Rhat = Rhat + fft(Gn_w_zeta(:,n-m+1,m+1)-Nn_zeta(:,n-m+1,m+1));
  end
  for m=0:n-1
    Qhat = Qhat - fft(Gn_u_U(:,n-m+1,m+1) + Gn_w_U(:,n-m+1,m+1)- Nn_U(:,n-m+1,m+1));
  end
  U_n(:,n+1) = ifft((Rhat + Qhat-fft(psi_n(:,n+1))).*Delta);

  % compute and store Gn_u[U_s] 
  xi = zeros(N_theta,N-s+1);
  xi(:,0+1) = U_n(:,s+1);
  anp = field_fe_helmholtz_polar_exterior(xi,f,k_u,a,p,N_theta,N-s);
  Gn = dno_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N-s);
  for r=0:N-s
    Gn_u_U(:,r+1,s+1) = Gn(:,r+1);
  end
  
  % compute and store Gn_w[U_s]
  xi = zeros(N_theta,N-s+1);
  xi(:,0+1) = U_n(:,s+1);
  dnp = field_fe_helmholtz_polar_interior(xi,f,k_w,a,p,N_theta,N-s);
  Gn = dno_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N-s);
  for r=0:N-s
    Gn_w_U(:,r+1,s+1) = Gn(:,r+1);
  end
  
  % compute and store |N|n[U_s]
  xi = zeros(N_theta,N-s+1);
  xi(:,0+1) = U_n(:,s+1);
  Nn = rho.*current_n(f_theta,xi,N-s);
  for r=0:N-s
    Nn_U(:,r+1,s+1) = Nn(:,r+1);
  end
end