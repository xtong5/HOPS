function [anp] = field_fe_helmholtz_upper_3d(xi_n,f,p1,p2,alpha1p,alpha2p,betap,...
    eep,eem,Nx1,Nx2,N)
anp = zeros(Nx1,Nx2,N+1); fn = zeros(Nx1,Nx2,N+1);
fn(:,:,0+1) = ones(Nx1,Nx2);
for n=1:N
  fn(:,:,n+1) = f.*fn(:,:,n-1+1)/n;
end
anp(:,:,0+1) = fft2( eem.*xi_n(:,:,0+1) );
for n=1:N
  anp(:,:,n+1) = fft2( eem.*xi_n(:,:,n+1) );
  for m=0:n-1
    anp(:,:,n+1) = anp(:,:,n+1) - fft2(eem.*( fn(:,:,n-m+1).*...
        eep.*ifft2((1i*betap).^(n-m).*anp(:,:,m+1)) ));
  end
end