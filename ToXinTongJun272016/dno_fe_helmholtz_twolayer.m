function [U_n] = dno_fe_helmholtz_twolayer(zeta_n,psi_n,f,tau2,...
    p,alphap,beta_up,beta_wp,eep,eem,Nx,N)
U_n = zeros(Nx,N+1);

Delta = -1i*beta_up - tau2*1i*beta_wp;
Rhat = -fft( eem.*psi_n(:,0+1) ) ...
    + tau2*(-1i*beta_wp).*fft( eem.*zeta_n(:,0+1) );
U_n(:,0+1) = eep.*ifft( Rhat./Delta );

for n=1:N
  Rhat = fft( -eem.*psi_n(:,n+1) );
  for m=0:n
    xi = zeros(Nx,n-m+1);
    xi(:,0+1) = zeta_n(:,m+1);
    dpn = field_fe_helmholtz_lower(xi,f,...
        p,alphap,beta_wp,eep,eem,Nx,n-m);
    Jn = dno_fe_helmholtz_lower(dpn,f,p,alphap,beta_wp,eep,eem,Nx,n-m);
    Rhat = Rhat + tau2*fft( eem.*Jn(:,n-m+1) );
  end
  
  for m=0:n-1
    xi = zeros(Nx,n-m+1);
    xi(:,0+1) = U_n(:,m+1);
    apn = field_fe_helmholtz_upper(xi,f,...
        p,alphap,beta_up,eep,eem,Nx,n-m);
    Gn = dno_fe_helmholtz_upper(apn,f,p,alphap,beta_up,eep,eem,Nx,n-m);
    Rhat = Rhat - fft( eem.*Gn(:,n-m+1) );
  end
  
  for m=0:n-1
    xi = zeros(Nx,n-m+1);
    xi(:,0+1) = U_n(:,m+1);
    dpn = field_fe_helmholtz_lower(xi,f,...
        p,alphap,beta_wp,eep,eem,Nx,n-m);
    Jn = dno_fe_helmholtz_lower(dpn,f,p,alphap,beta_wp,eep,eem,Nx,n-m);
    Rhat = Rhat - tau2*fft( eem.*Jn(:,n-m+1) );
  end
  
  U_n(:,n+1) = eep.*ifft( Rhat./Delta );
end