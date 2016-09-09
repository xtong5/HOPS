% n = 4;
% p = [0:8/2-1,-8/2:-1]';
% z = 2;
% a = diff_bessel(2,p,n,z);
% b = diff_besselh(p,n,z);
% a-b

L = 2*pi;
k_u = 3; 
k_w = 2;

Eps = 0.02;
N_theta = 64;
a = 2.0/5.0; 
N = 16;
theta = (L/N_theta)*[0:N_theta-1]';
p = [0:N_theta/2-1,-N_theta/2:-1]';
f = exp(cos(theta));
f_theta = -sin(theta).*f;
A=a+Eps.*f;

Ar_u = 2; r = 2; % take a special wavenumber
Ar_w = 1; %r = 2; % take a special wavenumber ???same

xi_u = Ar_u*besselh(r,k_u.*A).*exp(1i*r.*theta);
nu_u = Ar_u*((-k_u.*A.*(diff_bessel(2,r,1,k_u.*A))+...
    1i*r*Eps.*f_theta.*besselh(r,k_u.*A)./A ).*exp(1i*r.*theta)); 
xi_w = Ar_w*besselj(r,k_w.*A).*exp(1i*r.*theta);
nu_w = Ar_w*((k_w.*A.*(diff_bessel(1,r,1,k_w.*A))-...
    1i*r*Eps.*f_theta.*besselj(r,k_w.*A)./A ).*exp(1i*r.*theta)); 
xi_u_n = zeros(N_theta,N+1); nu_u_n = zeros(N_theta,N+1);
xi_w_n = zeros(N_theta,N+1); nu_w_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);f_nmt = ones(N_theta,1);
f_n = f.*f_n;
xi_u_n(:,0+1) = Ar_u*besselh(r,k_u*a).*exp(1i*r.*theta);
xi_u_n(:,1+1) = Ar_u*k_u^1*diff_bessel(2,r,1,k_u*a).*f_n.*exp(1i*r.*theta);
nu_u_n(:,0+1) = -Ar_u*k_u*a*diff_bessel(2,r,1,k_u*a).*exp(1i*r.*theta);
nu_u_n1(:,0+1) = -Ar_u*k_u*a*diff_besselh(r,1,k_u*a).*exp(1i*r.*theta);
nu_u_n(:,1+1) = -f/a.*nu_u_n(:,1)...
      -Ar_u*a*k_u^(1+1).*diff_bessel(2,r,1+1,k_u*a).*f_n.*exp(1i*r.*theta)...
      -Ar_u*(2*f).*k_u^1.*diff_bessel(2,r,1,k_u*a).*f_nmo.*exp(1i*r.*theta)...
      +Ar_u*(f_theta/a).*(1i*r).*besselh(r,k_u*a).*f_nmo.*exp(1i*r.*theta);
xi_w_n(:,0+1) = Ar_w*besselj(r,k_w*a).*exp(1i*r.*theta);
xi_w_n(:,1+1) = Ar_w*k_w^1*diff_bessel(1,r,1,k_w*a).*f_n.*exp(1i*r.*theta);
nu_w_n(:,0+1) = Ar_w*k_w*a*diff_bessel(1,r,1,k_w*a).*exp(1i*r.*theta);
nu_w_n1(:,0+1) = Ar_w*k_w*a*diff_besselj(r,1,k_w*a).*exp(1i*r.*theta);
nu_w_n(:,1+1) = -f/a.*nu_w_n(:,1)...
      +Ar_w*a*k_w^(1+1).*diff_bessel(1,r,1+1,k_w*a).*f_n.*exp(1i*r.*theta)...
      +Ar_w*(2*f).*k_w^1.*diff_bessel(1,r,1,k_w*a).*f_nmo.*exp(1i*r.*theta)...
      -Ar_w*(f_theta/a).*(1i*r).*besselj(r,k_w*a).*f_nmo.*exp(1i*r.*theta);

xi_u1 = Ar_u*besselh(r,k_u.*A).*exp(1i*r.*theta);
nu_u1 = Ar_u*((-k_u.*A.*(diff_besselh(r,1,k_u.*A))+...
    1i*r*Eps.*f_theta.*besselh(r,k_u.*A)./A ).*exp(1i*r.*theta)); 
xi_w1 = Ar_w*besselj(r,k_w.*A).*exp(1i*r.*theta);
nu_w1 = Ar_w*((k_w.*A.*(diff_besselj(r,1,k_w.*A))-...
    1i*r*Eps.*f_theta.*besselj(r,k_w.*A)./A ).*exp(1i*r.*theta)); 
xi_u_n1(:,0+1) = Ar_u*besselh(r,k_u*a).*exp(1i*r.*theta);
xi_u_n1(:,1+1) = Ar_u*k_u^1*diff_besselh(r,1,k_u*a).*f_n.*exp(1i*r.*theta);
nu_u_n1(:,1+1) = -f/a.*nu_u_n(:,1)...
      -Ar_u*a*k_u^(1+1).*diff_besselh(r,1+1,k_u*a).*f_n.*exp(1i*r.*theta)...
      -Ar_u*(2*f).*k_u^1.*diff_besselh(r,1,k_u*a).*f_nmo.*exp(1i*r.*theta)...
      +Ar_u*(f_theta/a).*(1i*r).*besselh(r,k_u*a).*f_nmo.*exp(1i*r.*theta);
xi_w_n1(:,0+1) = Ar_w*besselj(r,k_w*a).*exp(1i*r.*theta);
xi_w_n1(:,1+1) = Ar_w*k_w^1*diff_besselj(r,1,k_w*a).*f_n.*exp(1i*r.*theta);
nu_w_n1(:,1+1) = -f/a.*nu_w_n(:,1)...
      +Ar_w*a*k_w^(1+1).*diff_besselj(r,1+1,k_w*a).*f_n.*exp(1i*r.*theta)...
      +Ar_w*(2*f).*k_w^1.*diff_besselj(r,1,k_w*a).*f_nmo.*exp(1i*r.*theta)...
      -Ar_w*(f_theta/a).*(1i*r).*besselj(r,k_w*a).*f_nmo.*exp(1i*r.*theta);

  

for n=2:N
  f_n = f.*f_n/n;
  f_nmo = f.*f_nmo/(n-1);
  if(n>2)
    f_nmt = f.*f_nmt/(n-2);
  end
  xi_u_n(:,n+1) = Ar_u*k_u^n*diff_bessel(2,r,n,k_u*a).*f_n.*exp(1i*r.*theta);
  nu_u_n(:,n+1) = -f/a.*nu_u_n(:,n-1+1)...
      -Ar_u*a*k_u^(n+1).*diff_bessel(2,r,n+1,k_u*a).*f_n.*exp(1i*r.*theta)...
      -Ar_u*(2*f).*k_u^n.*diff_bessel(2,r,n,k_u*a).*f_nmo.*exp(1i*r.*theta)...
      -Ar_u*(f.^2/a)*k_u^(n-1).*diff_bessel(2,r,n-1,k_u*a).*f_nmt.*exp(1i*r.*theta)...
      +Ar_u*(f_theta/a)*k_u^(n-1).*(1i*r).*diff_bessel(2,r,n-1,k_u*a)...
      .*f_nmo.*exp(1i*r.*theta);
  xi_w_n(:,n+1) = Ar_w*k_w^n*diff_bessel(1,r,n,k_w*a).*f_n.*exp(1i*r.*theta);
  nu_w_n(:,n+1) = -f/a.*nu_w_n(:,n-1+1)...
      +Ar_w*a*k_w^(n+1).*diff_bessel(1,r,n+1,k_w*a).*f_n.*exp(1i*r.*theta)...
      +Ar_w*(2*f).*k_w^n.*diff_bessel(1,r,n,k_w*a).*f_nmo.*exp(1i*r.*theta)...
      +Ar_w*(f.^2/a)*k_w^(n-1).*diff_bessel(1,r,n-1,k_w*a).*f_nmt.*exp(1i*r.*theta)...
      -Ar_w*(f_theta/a)*k_w^(n-1).*(1i*r).*diff_bessel(1,r,n-1,k_w*a)...
      .*f_nmo.*exp(1i*r.*theta);
   xi_u_n1(:,n+1) = Ar_u*k_u^n*diff_besselh(r,n,k_u*a).*f_n.*exp(1i*r.*theta);
  nu_u_n1(:,n+1) = -f/a.*nu_u_n(:,n-1+1)...
      -Ar_u*a*k_u^(n+1).*diff_besselh(r,n+1,k_u*a).*f_n.*exp(1i*r.*theta)...
      -Ar_u*(2*f).*k_u^n.*diff_besselh(r,n,k_u*a).*f_nmo.*exp(1i*r.*theta)...
      -Ar_u*(f.^2/a)*k_u^(n-1).*diff_besselh(r,n-1,k_u*a).*f_nmt.*exp(1i*r.*theta)...
      +Ar_u*(f_theta/a)*k_u^(n-1).*(1i*r).*diff_besselh(r,n-1,k_u*a)...
      .*f_nmo.*exp(1i*r.*theta);
  xi_w_n1(:,n+1) = Ar_w*k_w^n*diff_besselj(r,n,k_w*a).*f_n.*exp(1i*r.*theta);
  nu_w_n1(:,n+1) = -f/a.*nu_w_n(:,n-1+1)...
      +Ar_w*a*k_w^(n+1).*diff_besselj(r,n+1,k_w*a).*f_n.*exp(1i*r.*theta)...
      +Ar_w*(2*f).*k_w^n.*diff_besselj(r,n,k_w*a).*f_nmo.*exp(1i*r.*theta)...
      +Ar_w*(f.^2/a)*k_w^(n-1).*diff_besselj(r,n-1,k_w*a).*f_nmt.*exp(1i*r.*theta)...
      -Ar_w*(f_theta/a)*k_w^(n-1).*(1i*r).*diff_besselj(r,n-1,k_w*a)...
      .*f_nmo.*exp(1i*r.*theta);
end




norm(nu_u-nu_u1);
norm(nu_w-nu_w1);
norm(nu_u_n-nu_u_n1)
norm(nu_w_n-nu_w_n1)
for i = 1:17
x(i) = norm(nu_u_n(:,i)-nu_u_n1(:,i),1);
y(i) = norm(nu_w_n(:,i)-nu_w_n1(:,i),1);
end





    