function [dnp_n] = field_fe_n_helmholtz_polar_interior(xi,fn,k,p,n,a,dnp)
dnp_n = 0*xi;
for m=0:n-1
    dnp_n = dnp_n - k^(n-m)*fft(fn(:,n-m+1).*...
        ifft(diff_besselj(p,n-m,k*a).*dnp(:,m+1)./besselj(p,k*a)));
end
end