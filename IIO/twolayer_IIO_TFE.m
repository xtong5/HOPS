function [I_u_n,I_w_n] = twolayer_IIO_TFE(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,b,c,N_theta,N,N_r,eta)

%% set up
I_u_n = zeros(N_theta,N+1);
I_w_n = zeros(N_theta,N+1);
Qn = zeros(N_theta,N+1,N+1); %IIO for Q_u_jI_u_i
Sn = zeros(N_theta,N+1,N+1); %IIO for S_w_jI_w_i
Q0_p = (-a*k_u * diff_besselh(p,1,k_u*a)/sigma_u - 1i*eta*besselh(p,k_u*a))./...
    (-a*k_u * diff_besselh(p,1,k_u*a)/sigma_u+1i*eta*besselh(p,k_u*a));
S0_p = (a*k_w * diff_besselj(p,1,k_w*a)/sigma_w - 1i*eta*besselj(p,k_w*a))./...
    (a*k_w * diff_besselj(p,1,k_w*a)/sigma_w + 1i*eta*besselj(p,k_w*a));
Delta = 1./(1-Q0_p.*S0_p);


%% compute order n=0
n = 0;
s = n;
Ahat = 1i*eta*fft( zeta_n(:,0+1) );
Bhat = fft( -psi_n(:,0+1) );
R1hat = Ahat + Bhat;
R2hat = -Ahat + Bhat;
I_u_n(:,0+1) = ifft((R1hat - S0_p.*R2hat).*Delta);
I_w_n(:,0+1) = ifft((-Q0_p.*R1hat + R2hat).*Delta);

% compute and store Qn[I_u_0]
xi_u = zeros(N_theta,N-s+1); xi_u(:,0+1) = I_u_n(:,s+1); 
[Un,Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_exterior(xi_u,f,f_theta,k_u,a,b,p,N_theta,N-s,N_r,sigma_u,eta);
Qn_tfe = IIO_tfe_helmholtz_polar_exterior(Un,Dr_Un,Dp_Un,f,f_theta,a,b,N-s,sigma_u,eta);
for r=0:N-s
  Qn(:,r+1,s+1) = Qn_tfe(:,r+1);
end

% compute and store Sn[I_w_0]
xi_w = zeros(N_theta,N-s+1); xi_w(:,0+1) = I_w_n(:,s+1);
[Wn,Dr_Wn,Dp_Wn] = field_tfe_IIO_helmholtz_polar_interior(xi_w,f,f_theta,k_w,a,c,p,N_theta,N-s,N_r,sigma_w,eta);
Sn_tfe = IIO_tfe_helmholtz_polar_interior(Wn,Dr_Wn,Dp_Wn,f,f_theta,a,c,N-s,sigma_w,eta);
for r=0:N-s
  Sn(:,r+1,s+1) = Sn_tfe(:,r+1);
end

%% order n>0
for n=1:N
  s=n;
  Ahat = 1i*eta*fft( zeta_n(:,n+1) );
  Bhat = fft( -psi_n(:,n+1) );
  R1hat = Ahat + Bhat;
  R2hat = -Ahat + Bhat;
  for m=0:n-1
    R1hat = R1hat - fft(Sn(:,n+1-m,m+1));
    R2hat = R2hat - fft(Qn(:,n+1-m,m+1));
  end
  I_u_n(:,n+1) = ifft((R1hat - S0_p.*R2hat).*Delta);
  I_w_n(:,n+1) = ifft((-Q0_p.*R1hat + R2hat).*Delta);

  % compute and store Qn[I_u_n]
  xi_u = zeros(N_theta,N-s+1); xi_u(:,0+1) = I_u_n(:,s+1); 
  [Un,Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_exterior(xi_u,f,f_theta,k_u,a,b,p,N_theta,N-s,N_r,sigma_u,eta);
  Qn_tfe = IIO_tfe_helmholtz_polar_exterior(Un,Dr_Un,Dp_Un,f,f_theta,a,b,N-s,sigma_u,eta);
  for r=0:N-s
    Qn(:,r+1,s+1) = Qn_tfe(:,r+1);
  end
  
  % compute and store Sn[I_w_n]  
  xi_w = zeros(N_theta,N-s+1); xi_w(:,0+1) = I_w_n(:,s+1);
  [Wn,Dr_Wn,Dp_Wn] = field_tfe_IIO_helmholtz_polar_interior(xi_w,f,f_theta,k_w,a,c,p,N_theta,N-s,N_r,sigma_w,eta);
  Sn_tfe = IIO_tfe_helmholtz_polar_interior(Wn,Dr_Wn,Dp_Wn,f,f_theta,a,c,N-s,sigma_w,eta);
  for r=0:N-s
    Sn(:,r+1,s+1) = Sn_tfe(:,r+1);
  end

end