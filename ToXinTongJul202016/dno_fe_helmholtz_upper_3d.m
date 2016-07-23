function [Gn] = dno_fe_helmholtz_upper_3d(anp,f,p1,p2,alpha1p,alpha2p,betap,...
    eep,eem,Nx1,Nx2,N)
Gn = zeros(Nx1,Nx2,N+1); fn = zeros(Nx1,Nx2,N+1);
fn(:,:,0+1) = ones(Nx1,Nx2);
for n=1:N
  fn(:,:,n+1) = f.*fn(:,:,n-1+1)/n;
end
f_x1 = ifft2( (1i*p1).*fft2(f) );
f_x2 = ifft2( (1i*p2).*fft2(f) );
for n=0:N
  Gn(:,:,n+1) = eep.*ifft2( -(1i*betap).*anp(:,:,n+1) );
  for m=0:n-1
    Gn(:,:,n+1) = Gn(:,:,n+1) - fn(:,:,n-m+1).*...
        eep.*ifft2( (1i*betap).^(n+1-m).*anp(:,:,m+1) )...
        + f_x1.*fn(:,:,n-1-m+1).*eep.*ifft2( (1i*alpha1p).*...
        (1i*betap).^(n-1-m).*anp(:,:,m+1) )...
        + f_x2.*fn(:,:,n-1-m+1).*eep.*ifft2( (1i*alpha2p).*...
        (1i*betap).^(n-1-m).*anp(:,:,m+1) );
  end
end