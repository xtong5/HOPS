function [anp,dnp] = field_fe_helmholtz_twolayer(zeta_n,psi_n,f,tau2,...
    p,alphap,beta_up,beta_wp,eep,eem,Nx,N)
anp = zeros(Nx,N+1); dnp = zeros(Nx,N+1); fn = zeros(Nx,N+1);
fn(:,0+1) = ones(Nx,1);
for n=1:N
  fn(:,n+1) = f.*fn(:,n-1+1)/n;
end
f_x = real(ifft( (1i*p).*fft(f) ));
Delta = tau2*1i*beta_wp + 1i*beta_up;
M0Inv_11 = tau2*(1i*beta_wp)./Delta;
M0Inv_12 = 1.0./Delta;
M0Inv_21 = -1i*beta_up./Delta;
M0Inv_22 = 1.0./Delta;
Qhat = fft( eem.*zeta_n(:,0+1) );
Rhat = fft( eem.*psi_n(:,0+1) );
anp(:,0+1) = M0Inv_11.*Qhat + M0Inv_12.*Rhat;
dnp(:,0+1) = M0Inv_21.*Qhat + M0Inv_22.*Rhat;

for n=1:N
  Qhat = fft( eem.*zeta_n(:,n+1) );
  Rhat = fft( eem.*psi_n(:,n+1) );
  for m=0:n-1
    Qhat = Qhat - fft(eem.*( fn(:,n-m+1).*...
        eep.*ifft( (1i*beta_up).^(n-m).*anp(:,m+1)) ))...
        + fft(eem.*( fn(:,n-m+1).*...
        eep.*ifft( (-1i*beta_wp).^(n-m).*dnp(:,m+1)) ));
    Rhat = Rhat - fft(eem.*( fn(:,n-m+1).*...
        eep.*ifft( (1i*beta_up).^(n-m+1).*anp(:,m+1)) ))...
        + tau2*fft(eem.*( fn(:,n-m+1).*...
        eep.*ifft( (-1i*beta_wp).^(n-m+1).*dnp(:,m+1)) ));

    Rhat = Rhat + fft(eem.*( f_x.*fn(:,n-m-1+1).*...
      eep.*ifft( (1i*alphap).*(1i*beta_up).^(n-m-1).*anp(:,m+1)) ))...
      - tau2*fft(eem.*( f_x.*fn(:,n-m-1+1).*...
      eep.*ifft( (1i*alphap).*(-1i*beta_wp).^(n-m-1).*dnp(:,m+1)) ));
  end
  anp(:,n+1) = M0Inv_11.*Qhat + M0Inv_12.*Rhat;
  dnp(:,n+1) = M0Inv_21.*Qhat + M0Inv_22.*Rhat;
end