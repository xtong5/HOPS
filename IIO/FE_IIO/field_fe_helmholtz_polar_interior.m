function [dnp] = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k,a,p,N_theta,N,sigma,Z_p)
dnp = zeros(N_theta,N+1); 
fn = zeros(N_theta,N+1);
fn(:,0+1) = ones(N_theta,1);
Delta = 1./(sigma*a*k*diff_besselj(p,1,k*a)./besselj(p,k*a) - Z_p);
for n=1:N
  fn(:,n+1) = f.*fn(:,n-1+1)/n; %Fn
end
dnp(:,0+1) = fft(I_w_n(:,0+1)).*Delta; 

for n=1:N
    Ln = I_w_n(:,n+1);
    for m=0:n-1
      Ln = Ln - sigma*(a.*ifft(k^(n-m+1).*diff_besselj(p,n-m+1,k*a).*dnp(:,m+1)./besselj(p,k*a)).*fn(:,n-m+1)...
          +2*f.*ifft(k^(n-m).*diff_besselj(p,n-m,k*a).*dnp(:,m+1)./besselj(p,k*a)).*fn(:,n-1-m+1)...
          -f_theta./a.* ifft(1i.*p.*k^(n-m-1).*diff_besselj(p,n-m-1,k*a).*dnp(:,m+1)./besselj(p,k*a)).*fn(:,n-1-m+1))...
          + ifft(Z_p.*fft(ifft(k^(n-m).*diff_besselj(p,n-m,k*a).*dnp(:,m+1)./besselj(p,k*a)).*fn(:,n-m+1)))...
          + f./a.*(ifft(Z_p.*fft(ifft(k^(n-m-1).*diff_besselj(p,n-m-1,k*a).*dnp(:,m+1)./besselj(p,k*a)).*fn(:,n-1-m+1))));
    end
    for m=0:n-2
      Ln = Ln - sigma*f.^2.*fn(:,n-2-m+1).*ifft(k^(n-m-1).*diff_besselj(p,n-m-1,k*a).*dnp(:,m+1)./besselj(p,k*a))/a;
    end
    Ln = Ln + f./a.*I_w_n(:,n-1+1);
    dnp(:,n+1) = fft(Ln).*Delta;
end
end