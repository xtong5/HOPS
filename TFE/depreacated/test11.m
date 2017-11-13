Ar = 2; pp = 2; % compute a special wavenumber
AA=a+Eps.*f;
xi = Ar*besselh(pp,k.*AA).*exp(1i*pp.*theta);
nu1 =Ar*((-k.*AA.*(diff_bessel(2,pp,1,k.*AA))+...
    1i*pp*Eps.*f_theta.*besselh(pp,k.*AA)./AA ).*exp(1i*pp.*theta)); %DNO
xi_n = zeros(N_theta,N+1); nu_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1); f_nmt = ones(N_theta,1);
xi_n(:,0+1) = Ar*besselh(pp,k*a).*exp(1i*pp.*theta);
nu_n(:,0+1) = -Ar*k*a*diff_bessel(2,pp,1,k*a).*exp(1i*pp.*theta);
f_n = f.*f_n;
xi_n(:,1+1) = Ar*k^1*diff_bessel(2,pp,1,k*a).*f_n.*exp(1i*pp.*theta);
nu_n(:,1+1) = -f/a.*nu_n(:,1)...
      -Ar*a*k^(1+1).*diff_bessel(2,pp,1+1,k*a).*f_n.*exp(1i*pp.*theta)...
      -Ar*(2*f).*k^1.*diff_bessel(2,pp,1,k*a).*f_nmo.*exp(1i*pp.*theta)...
      +Ar*(f_theta/a).*(1i*pp).*besselh(pp,k*a)...
      .*f_nmo.*exp(1i*pp.*theta);

for n=2:N
  f_n = f.*f_n/n;
  f_nmo = f.*f_nmo/(n-1);

  if(n>2)
    f_nmt = f.*f_nmt/(n-2);
  end
  xi_n(:,n+1) = Ar*k^n*diff_bessel(2,pp,n,k*a).*f_n.*exp(1i*pp.*theta);
  nu_n(:,n+1) = -f/a.*nu_n(:,n-1+1)...
      -Ar*a*k^(n+1).*diff_bessel(2,pp,n+1,k*a).*f_n.*exp(1i*pp.*theta)...
      -Ar*(2*f).*k^n.*diff_bessel(2,pp,n,k*a).*f_nmo.*exp(1i*pp.*theta)...
      -Ar*(f.^2/a)*k^(n-1).*diff_bessel(2,pp,n-1,k*a).*f_nmt.*exp(1i*pp.*theta)...
      +Ar*(f_theta/a)*k^(n-1).*(1i*pp).*diff_bessel(2,pp,n-1,k*a)...
      .*f_nmo.*exp(1i*pp.*theta);
end

nn=16;
norm(nu_n(:,nn)-Gn_tfe(:,nn))
     

