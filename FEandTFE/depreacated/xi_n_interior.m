function [xi_n_p] = xi_n_interior(N_theta,N,f,Ar,pp,a,k,theta)
f_n = ones(N_theta,1); 
xi_n_p(:,0+1) = fft(Ar*besselj(pp,k*a).*exp(1i*pp.*theta));
for n=1:N
  f_n = f.*f_n/n;
  xi_n_p(:,n+1) = fft(Ar*k^n*diff_besselj(pp,n,k*a).*f_n.*exp(1i*pp.*theta));
end

end