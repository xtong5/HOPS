function [anp_n] = field_fe_n_helmholtz_polar_exterior(xi,fn,k,p,n,a,anp)
anp_n = 0*xi;
for m=0:n-1
    anp_n = anp_n - k^(n-m)*fft(fn(:,n-m+1).*...
        ifft(diff_besselh(p,n-m,k*a).*anp(:,m+1)./besselh(p,k*a)));
end
end