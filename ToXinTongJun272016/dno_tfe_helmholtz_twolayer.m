function [U_n] = dno_tfe_helmholtz_twolayer(zeta_n,psi_n,f,tau2,...
    p,alphap,beta_up,beta_wp,eep,eem,Dy,a,b,Nx,Ny,N)
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
    wn = field_tfe_helmholtz_lower(xi,f,...
        p,alphap,beta_wp,eep,eem,Dy,b,Nx,Ny,n-m);
    Jn = dno_tfe_helmholtz_lower(wn,f,...
        p,alphap,beta_wp,eep,eem,Dy,b,Nx,Ny,n-m);
    Rhat = Rhat + tau2*fft( eem.*Jn(:,n-m+1) );
  end
  
  for m=0:n-1
    xi = zeros(Nx,n-m+1);
    xi(:,0+1) = U_n(:,m+1);
    un = field_tfe_helmholtz_upper(xi,f,...
        p,alphap,beta_up,eep,eem,Dy,a,Nx,Ny,n-m);
    Gn = dno_tfe_helmholtz_upper(un,f,...
        p,alphap,beta_up,eep,eem,Dy,a,Nx,Ny,n-m);
    Rhat = Rhat - fft( eem.*Gn(:,n-m+1) );
  end
  
  for m=0:n-1
    xi = zeros(Nx,n-m+1);
    xi(:,0+1) = U_n(:,m+1);
    wn = field_tfe_helmholtz_lower(xi,f,...
        p,alphap,beta_wp,eep,eem,Dy,b,Nx,Ny,n-m);
    Jn = dno_tfe_helmholtz_lower(wn,f,...
        p,alphap,beta_wp,eep,eem,Dy,b,Nx,Ny,n-m);
    Rhat = Rhat - tau2*fft( eem.*Jn(:,n-m+1) );
  end
  
  U_n(:,n+1) = eep.*ifft( Rhat./Delta );
end