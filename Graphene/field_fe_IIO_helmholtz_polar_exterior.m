function [anp] = field_fe_IIO_helmholtz_polar_exterior(I_u_n,f,f_theta,k,a,p,N_theta,N,sigma,Y_p)
%check the order of input
anp = zeros(N_theta,N+1); 
fn = zeros(N_theta,N+1);
fn(:,0+1) = ones(N_theta,1);
Delta = 1./(-sigma*a*k*diff_besselh(p,1,k*a)./besselh(p,k*a)+ Y_p);
for n=1:N
  fn(:,n+1) = f.*fn(:,n-1+1)/n; %Fn
end
anp(:,0+1) = fft(I_u_n(:,0+1)).*Delta; 

for n=1:N
    Ln = I_u_n(:,n+1);
    for m=0:n-1
      Ln = Ln + sigma*(a.*ifft(k^(n-m+1).*diff_besselh(p,n-m+1,k*a).*anp(:,m+1)./besselh(p,k*a)).*fn(:,n-m+1)...
          +2*f.*ifft(k^(n-m).*diff_besselh(p,n-m,k*a).*anp(:,m+1)./besselh(p,k*a)).*fn(:,n-1-m+1)...
          -f_theta./a.* ifft(1i.*p.*k^(n-m-1).*diff_besselh(p,n-m-1,k*a).*anp(:,m+1)./besselh(p,k*a)).*fn(:,n-1-m+1))...
          -ifft(Y_p.*fft(ifft(k^(n-m).*diff_besselh(p,n-m,k*a).*anp(:,m+1)./besselh(p,k*a)).*fn(:,n-m+1)))...
          -f./a.*(ifft(Y_p.*fft(ifft(k^(n-m-1).*diff_besselh(p,n-m-1,k*a).*anp(:,m+1)./besselh(p,k*a)).*fn(:,n-1-m+1))));
    end
    for m=0:n-2
      Ln = Ln + sigma*f.^2.*fn(:,n-m-1).*ifft(k^(n-m-1).*diff_besselh(p,n-m-1,k*a).*anp(:,m+1)./besselh(p,k*a))/a;
    end
    Ln = Ln + f./a.*I_u_n(:,n-1+1);
    anp(:,n+1) = fft(Ln).*Delta;
end
end