function [U_n] = twolayer_dno_fe_helmholtz(zeta_n,psi_n,f,tau2,...
    p,alphap,beta_up,beta_wp,eep,eem,Nx,N)
U_n = zeros(Nx,N+1);
J_zeta_rs = zeros(Nx,N+1,N+1);
G_U_rs = zeros(Nx,N+1,N+1);
J_U_rs = zeros(Nx,N+1,N+1);

% Compute and store J_r[zeta_s]
for n=0:N
  s = n;
  xi = zeros(Nx,N-s+1);
  xi(:,0+1) = zeta_n(:,s+1);
  dpn = field_fe_helmholtz_lower(xi,f,p,alphap,beta_wp,eep,eem,Nx,N-s);
  Jn = dno_fe_helmholtz_lower(dpn,f,p,alphap,beta_wp,eep,eem,Nx,N-s);
  for r=0:N-s
    J_zeta_rs(:,r+1,s+1) = Jn(:,r+1);
  end
end

%
% n = 0
%

n = 0;
Delta = -1i*beta_up - tau2*1i*beta_wp;
Rhat = -fft( eem.*psi_n(:,0+1) ) ...
    + tau2*(-1i*beta_wp).*fft( eem.*zeta_n(:,0+1) );
U_n(:,0+1) = eep.*ifft( Rhat./Delta );

% Compute and store G_r[U_s]
s = n;
xi = zeros(Nx,N-s+1);
xi(:,0+1) = U_n(:,s+1);
apn = field_fe_helmholtz_upper(xi,f,p,alphap,beta_up,eep,eem,Nx,N-s);
Gn = dno_fe_helmholtz_upper(apn,f,p,alphap,beta_up,eep,eem,Nx,N-s);
for r=0:N-s
  G_U_rs(:,r+1,s+1) = Gn(:,r+1);
end

% Compute and store J_r[U_s]
s = n;
xi = zeros(Nx,N-s+1);
xi(:,0+1) = U_n(:,s+1);
dpn = field_fe_helmholtz_lower(xi,f,p,alphap,beta_wp,eep,eem,Nx,N-s);
Jn = dno_fe_helmholtz_lower(dpn,f,p,alphap,beta_wp,eep,eem,Nx,N-s);
for r=0:N-s
  J_U_rs(:,r+1,s+1) = Jn(:,r+1);
end

%
% Order n>0
%

for n=1:N
  Rhat = fft( -eem.*psi_n(:,n+1) );
  for m=0:n
     Rhat = Rhat + tau2*fft( eem.*J_zeta_rs(:,n-m+1,m+1) );
   end
  
  for m=0:n-1
     Rhat = Rhat - fft( eem.*G_U_rs(:,n-m+1,m+1) );
   end
  
  for m=0:n-1
     Rhat = Rhat - tau2*fft( eem.*J_U_rs(:,n-m+1,m+1) );
   end
  
  U_n(:,n+1) = eep.*ifft( Rhat./Delta );
  
  % Compute and store G_r[U_s]
  s = n;
  xi = zeros(Nx,N-s+1);
  xi(:,0+1) = U_n(:,s+1);
  apn = field_fe_helmholtz_upper(xi,f,p,alphap,beta_up,eep,eem,Nx,N-s);
  Gn = dno_fe_helmholtz_upper(apn,f,p,alphap,beta_up,eep,eem,Nx,N-s);
  for r=0:N-s
    G_U_rs(:,r+1,s+1) = Gn(:,r+1);
  end

  % Compute and store J_r[U_s]
  s = n;
  xi = zeros(Nx,N-s+1);
  xi(:,0+1) = U_n(:,s+1);
  dpn = field_fe_helmholtz_lower(xi,f,p,alphap,beta_wp,eep,eem,Nx,N-s);
  Jn = dno_fe_helmholtz_lower(dpn,f,p,alphap,beta_wp,eep,eem,Nx,N-s);
  for r=0:N-s
    J_U_rs(:,r+1,s+1) = Jn(:,r+1);
  end
end