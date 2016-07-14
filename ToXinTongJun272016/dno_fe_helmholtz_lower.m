function [Gn] = dno_fe_helmholtz_lower(anp,f,p,alphap,betap,eep,eem,Nx,N)
Gn = zeros(Nx,N+1); fn = zeros(Nx,N+1);
fn(:,0+1) = ones(Nx,1);
for n=1:N
  fn(:,n+1) = f.*fn(:,n-1+1)/n;
end
f_x = ifft( (1i*p).*fft(f) );
for n=0:N
  Gn(:,n+1) = eep.*ifft( (-1i*betap).*anp(:,n+1) );
  for m=0:n-1
    Gn(:,n+1) = Gn(:,n+1) + fn(:,n-m+1).*...
        eep.*ifft( (-1i*betap).^(n+1-m).*anp(:,m+1) )...
        - f_x.*fn(:,n-1-m+1).*eep.*ifft( (1i*alphap).*...
        (-1i*betap).^(n-1-m).*anp(:,m+1) );
  end
end