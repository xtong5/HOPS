function [dnp] = field_fe_helmholtz_polar(xi,f,k,p,N_theta,N,a)
dnp = zeros(N_theta,N+1); 
fn = zeros(N_theta,N+1);
fn(:,0+1) = ones(N_theta,1);
for n=1:N
  fn(:,n+1) = f.*fn(:,n-1+1)/n; %Fn
end
dnp(:,0+1) = fft(xi); 
for n=1:N
  dnp(:,n+1) = 0*xi;
  for m=0:n-1
    dnp(:,n+1) = dnp(:,n+1) - k^(n-m+1)*fft(fn(:,n-m+1).*...
        ifft(diff_besselh(p,n-m+1,k*a).*dnp(:,m+1)./besselh(p,k*a))); 
  end
end