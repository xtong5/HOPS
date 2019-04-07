function [anp] = field_fe_helmholtz_polar_exterior(xi_n,f,k,a,p,N_theta,N)
%check the order of input
anp = zeros(N_theta,N+1); 
fn = zeros(N_theta,N+1);
fn(:,0+1) = ones(N_theta,1);
for n=1:N
  fn(:,n+1) = f.*fn(:,n-1+1)/n; %Fn
end
anp(:,0+1) = fft(xi_n(:,0+1)); 
for n=1:N
  anp(:,n+1) = fft(xi_n(:,n+1)); 
  for m=0:n-1
    anp(:,n+1) = anp(:,n+1) - k^(n-m)*fft(fn(:,n-m+1).*...
        ifft(diff_bessel(2,p,n-m,k*a).*anp(:,m+1)./besselh(p,k*a)));
  end
end