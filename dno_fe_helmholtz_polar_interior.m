function [Gn] = dno_fe_helmholtz_polar_interior(dnp,f,f_theta,k,a,p,N_theta,N)
Gn = zeros(N_theta,N+1); fn = zeros(N_theta,N+1);
fn(:,0+1) = ones(N_theta,1);
z=k*a;
for n=1:N
  fn(:,n+1) = f.*fn(:,n-1+1)/n;
end
%f_theta = ifft( (1i*p).*fft(f) ); 
Gn(:,1) = a*ifft(k.*diff_besselj(p,1,z).*dnp(:,1)./(besselj(p.',z).'));
for n=1:N
    Gn(:,n+1) = -f.* Gn(:,n)/a;
    for m=0:n
        Gn(:,n+1) = Gn(:,n+1)+a*fn(:,n-m+1).*ifft(k^(n-m+1).*...
            diff_besselj(p,n-m+1,z).*dnp(:,m+1)./(besselj(p.',z).'));
    end
    for m=0:n-2
        Gn(:,n+1) = Gn(:,n+1)+f.^2.*fn(:,n-m-1).*ifft(k^(n-m-1).*...
            diff_besselj(p,n-m-1,z).*dnp(:,m+1)./(besselj(p.',z).'))/a;
    end
    for m=0:n-1
        Gn(:,n+1) = Gn(:,n+1)+2*f.*fn(:,n-m).*ifft(k^(n-m).*...
            diff_besselj(p,n-m,z).*dnp(:,m+1)./(besselj(p.',z).'))...
            -f_theta.*fn(:,n-m).*ifft(k^(n-m).*(1i*p).*...
            diff_besselj(p,n-m-1,z).*dnp(:,m+1)./(besselj(p.',z).'))/a;
    end
end
end