function [Qn] = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k,a,p,N_theta,N,sigma,Z_p)
Qn = zeros(N_theta,N+1); fn = zeros(N_theta,N+1);
fn(:,0+1) = ones(N_theta,1);
z=k*a;
for n=1:N
  fn(:,n+1) = f.*fn(:,n-1+1)/n;
end
Qn(:,1) = -sigma*a*ifft(k.*diff_besselh(p,1,z).*anp(:,1)./besselh(p,z))+ifft(Z_p.*anp(:,1));
for n=1:N
    Qn(:,n+1) = -f.* Qn(:,n)/a;
    for m=0:n
        Qn(:,n+1) = Qn(:,n+1)-sigma*a*fn(:,n-m+1).*ifft(k^(n-m+1).*...
            diff_besselh(p,n-m+1,z).*anp(:,m+1)./besselh(p,z))...
            + ifft(Z_p.*fft(fn(:,n-m+1).*ifft(k^(n-m).*diff_besselh(p,n-m,z).*anp(:,m+1)./besselh(p,z))));
    end
    for m=0:n-2
        Qn(:,n+1) = Qn(:,n+1)-sigma*f.^2.*fn(:,n-m-1).*ifft(k^(n-m-1).*...
            diff_besselh(p,n-m-1,z).*anp(:,m+1)./besselh(p,z))/a;
    end
    for m=0:n-1
        Qn(:,n+1) = Qn(:,n+1)-sigma*2*f.*fn(:,n-m).*ifft(k^(n-m).*...
            diff_besselh(p,n-m,z).*anp(:,m+1)./besselh(p,z))...
            +sigma*f_theta.*fn(:,n-m).*ifft(k^(n-m-1).*(1i*p).*...
            diff_besselh(p,n-m-1,z).*anp(:,m+1)./besselh(p,z))/a...
            +f.*ifft(Z_p.*fft(fn(:,n-m).*ifft(k^(n-m-1).*diff_besselh(p,n-m-1,z).*anp(:,m+1)./besselh(p,z))))/a;
    end
end
end