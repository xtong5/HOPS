function [U_n] = twolayer_dno_tfe_helmholtz_polar(zeta_n,psi_n,f,f_theta,tau2,...
    p,k_u,k_w,a,b,c,N_theta,N,N_r)

%% set up
U_n = zeros(N_theta,N+1);
Gn_u_U = zeros(N_theta,N+1,N+1); %dno for G_u_jU_i
Gn_w_U = zeros(N_theta,N+1,N+1); %dno for G_w_jU_i
Gn_w_zeta = zeros(N_theta,N+1,N+1); %dno for G_w_jzeta_i
Delta1 = -a*k_u*diff_bessel(2,p,1,a*k_u)./besselh(p,a*k_u);
Delta2 = a*k_w*diff_bessel(1,p,1,a*k_w)./besselj(p,a*k_w);
Delta = Delta1 + tau2*Delta2;


%% Compute and store Gn_w[zeta]
for n=0:N
    s = n;
    xi = zeros(N_theta,N-s+1);
    xi(:,0+1) = zeta_n(:,s+1);
    [Un,Dr_Un,Dp_Un] = field_tfe_helmholtz_polar_interior(xi,f,f_theta,k_w,a,c,p,N_theta,N-s,N_r);
    Gn = dno_tfe_helmholtz_polar_interior(Dr_Un,Dp_Un,f,f_theta,k_w,a,c,p,N_theta,N-s,N_r);
    for r=0:N-s
        Gn_w_zeta(:,r+1,s+1) = Gn(:,r+1);
    end
end


%% compute order n=0
n = 0;
s = n;
Qhat = fft( zeta_n(:,0+1) );
Rhat = fft( -psi_n(:,0+1) );
U_n(:,0+1) = ifft((Rhat + tau2.*Delta2.*Qhat)./Delta);

% compute and store Gn_u[U_0]
xi = zeros(N_theta,N-s+1);
xi(:,0+1) = U_n(:,s+1);
[Un,Dr_Un,Dp_Un] = field_tfe_helmholtz_polar_exterior(xi,f,f_theta,k_u,a,b,p,N_theta,N-s,N_r);
Gn = dno_tfe_helmholtz_polar_exterior(Dr_Un,Dp_Un,f,f_theta,k_u,a,b,p,N_theta,N-s,N_r);
for r=0:N-s
  Gn_u_U(:,r+1,s+1) = Gn(:,r+1);
end

% compute and store Gn_w[U_0]
xi = zeros(N_theta,N-s+1);
xi(:,0+1) = U_n(:,s+1);
[Un,Dr_Un,Dp_Un] = field_tfe_helmholtz_polar_interior(xi,f,f_theta,k_w,a,c,p,N_theta,N-s,N_r);
Gn = dno_tfe_helmholtz_polar_interior(Dr_Un,Dp_Un,f,f_theta,k_w,a,c,p,N_theta,N-s,N_r);
for r=0:N-s
  Gn_w_U(:,r+1,s+1) = Gn(:,r+1);
end

%% order n>0
for n=1:N
  s=n;
  Rhat = fft(-psi_n(:,n+1) );
  Qhat = 0;
  for m=0:n
    Rhat = Rhat + tau2*fft(Gn_w_zeta(:,n+1-m,m+1));
  end
  for m=0:n-1
    Qhat = Qhat -fft(Gn_u_U(:,n+1-m,m+1)) - tau2*fft(Gn_w_U(:,n+1-m,m+1));
  end
  U_n(:,n+1) = ifft((Rhat + Qhat)./Delta);

  % compute and store Gn_u[U_s] 
  xi = zeros(N_theta,N-s+1);
  xi(:,0+1) = U_n(:,s+1);
  [Un,Dr_Un,Dp_Un] = field_tfe_helmholtz_polar_exterior(xi,f,f_theta,k_u,a,b,p,N_theta,N-s,N_r);
  Gn = dno_tfe_helmholtz_polar_exterior(Dr_Un,Dp_Un,f,f_theta,k_u,a,b,p,N_theta,N-s,N_r);
  for r=0:N-s
    Gn_u_U(:,r+1,s+1) = Gn(:,r+1);
  end
  
  % compute and store Gn_w[U_s]
  xi = zeros(N_theta,N-s+1);
  xi(:,0+1) = U_n(:,s+1);
  [Un,Dr_Un,Dp_Un] = field_tfe_helmholtz_polar_interior(xi,f,f_theta,k_w,a,c,p,N_theta,N-s,N_r);
  Gn = dno_tfe_helmholtz_polar_interior(Dr_Un,Dp_Un,f,f_theta,k_w,a,c,p,N_theta,N-s,N_r);
  for r=0:N-s
    Gn_w_U(:,r+1,s+1) = Gn(:,r+1);
  end
end