% check the GNO
Mode = 1;

L = 2*pi;
k_u = 1; 
k_w = 1;

Eps = 0.02;
N_theta = 64;
a = 1; 
N = 16;

fprintf('test_helmholtz_twolayer_polar\n');
fprintf('-------------\n');
fprintf('RunNumber = %d\n',RunNumber);
fprintf('k_u = %g  k_w = %g\n\n',k_u,k_w);
fprintf('Eps = %g  a = %g\n',Eps,a);
fprintf('N_theta = %d N = %d\n',N_theta,N);
fprintf('\n');

theta = (L/N_theta)*[0:N_theta-1]';
p = [0:N_theta/2-1,-N_theta/2:-1]';

f = exp(cos(theta));
f_theta = -sin(theta).*f;
A=a+Eps.*f;

Ar_u = 2; r = 2; % take a special wavenumber
Ar_w = 1; %r = 2; % take a special wavenumber ???same

xi_u = Ar_u*besselh(r,k_u.*A).*exp(1i*r.*theta);
nu_u = Ar_u*((-k_u.*A.*(diff_besselh(r,1,k_u.*A))+...
    1i*r*Eps.*f_theta.*besselh(r,k_u.*A)./A ).*exp(1i*r.*theta)); 
xi_w = Ar_w*besselj(r,k_w.*A).*exp(1i*r.*theta);
nu_w = Ar_w*((k_w.*A.*(diff_besselj(r,1,k_w.*A))-...
    1i*r*Eps.*f_theta.*besselj(r,k_w.*A)./A ).*exp(1i*r.*theta)); 
xi_u_n = zeros(N_theta,N+1); nu_u_n = zeros(N_theta,N+1);
xi_w_n = zeros(N_theta,N+1); nu_w_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);f_nmt = ones(N_theta,1);
f_n = f.*f_n;
xi_u_n(:,0+1) = Ar_u*besselh(r,k_u*a).*exp(1i*r.*theta);
xi_u_n(:,1+1) = Ar_u*k_u^1*diff_besselh(r,1,k_u*a).*f_n.*exp(1i*r.*theta);
nu_u_n(:,0+1) = -Ar_u*k_u*a*diff_besselh(r,1,k_u*a).*exp(1i*r.*theta);
nu_u_n(:,1+1) = -f/a.*nu_u_n(:,1)...
      -Ar_u*a*k_u^(1+1).*diff_besselh(r,1+1,k_u*a).*f_n.*exp(1i*r.*theta)...
      -Ar_u*(2*f).*k_u^1.*diff_besselh(r,1,k_u*a).*f_nmo.*exp(1i*r.*theta)...
      +Ar_u*(f_theta/a).*(1i*r).*besselh(r,k_u*a).*f_nmo.*exp(1i*r.*theta);
xi_w_n(:,0+1) = Ar_w*besselj(r,k_w*a).*exp(1i*r.*theta);
xi_w_n(:,1+1) = Ar_w*k_w^1*diff_besselj(r,1,k_w*a).*f_n.*exp(1i*r.*theta);
nu_w_n(:,0+1) = Ar_w*k_w*a*diff_besselj(r,1,k_w*a).*exp(1i*r.*theta);
nu_w_n(:,1+1) = -f/a.*nu_w_n(:,1)...
      +Ar_w*a*k_w^(1+1).*diff_besselj(r,1+1,k_w*a).*f_n.*exp(1i*r.*theta)...
      +Ar_w*(2*f).*k_w^1.*diff_besselj(r,1,k_w*a).*f_nmo.*exp(1i*r.*theta)...
      -Ar_w*(f_theta/a).*(1i*r).*besselj(r,k_w*a).*f_nmo.*exp(1i*r.*theta);

for n=2:N
  f_n = f.*f_n/n;
  f_nmo = f.*f_nmo/(n-1);
  if(n>2)
    f_nmt = f.*f_nmt/(n-2);
  end
  xi_u_n(:,n+1) = Ar_u*k_u^n*diff_besselh(r,n,k_u*a).*f_n.*exp(1i*r.*theta);
  nu_u_n(:,n+1) = -f/a.*nu_u_n(:,n-1+1)...
      -Ar_u*a*k_u^(n+1).*diff_besselh(r,n+1,k_u*a).*f_n.*exp(1i*r.*theta)...
      -Ar_u*(2*f).*k_u^n.*diff_besselh(r,n,k_u*a).*f_nmo.*exp(1i*r.*theta)...
      -Ar_u*(f.^2/a)*k_u^(n-1).*diff_besselh(r,n-1,k_u*a).*f_nmt.*exp(1i*r.*theta)...
      +Ar_u*(f_theta/a)*k_u^(n-1).*(1i*r).*diff_besselh(r,n-1,k_u*a)...
      .*f_nmo.*exp(1i*r.*theta);
  xi_w_n(:,n+1) = Ar_w*k_w^n*diff_besselj(r,n,k_w*a).*f_n.*exp(1i*r.*theta);
  nu_w_n(:,n+1) = -f/a.*nu_w_n(:,n-1+1)...
      +Ar_w*a*k_w^(n+1).*diff_besselj(r,n+1,k_w*a).*f_n.*exp(1i*r.*theta)...
      +Ar_w*(2*f).*k_w^n.*diff_besselj(r,n,k_w*a).*f_nmo.*exp(1i*r.*theta)...
      +Ar_w*(f.^2/a)*k_w^(n-1).*diff_besselj(r,n-1,k_w*a).*f_nmt.*exp(1i*r.*theta)...
      -Ar_w*(f_theta/a)*k_w^(n-1).*(1i*r).*diff_besselj(r,n-1,k_w*a)...
      .*f_nmo.*exp(1i*r.*theta);
end

% testing
%xi_u_n = zeros(N_theta,N+1);
%nu_u_n = zeros(N_theta,N+1);
%xi_w_n = zeros(N_theta,N+1);
%nu_w_n = zeros(N_theta,N+1);
%xi_u_n(:,0+1) = xi_u;
%nu_u_n(:,0+1) = nu_u;
%xi_w_n(:,0+1) = xi_w;
%nu_w_n(:,0+1) = nu_w;
for n=1:N
  for j=1:N_theta
    nu_u_taylor(j) = taylorsum(nu_u_n(j,:)',Eps,n);
    nu_w_taylor(j) = taylorsum(nu_w_n(j,:)',Eps,n);
  end
end 

nu_u_sum = zeros(N_theta,1);
nu_w_sum = zeros(N_theta,1);
for n=1:N+1
    EPS = Eps*(n-1);
    nu_u_sum = nu_u_sum + nu_u_n(:,n).* EPS;
    nu_w_sum = nu_w_sum + nu_w_n(:,n).* EPS;
end


nu_u_app = norm(nu_u-nu_u_taylor')
nu_w_app = norm(nu_w-nu_w_taylor')

error1 = norm(nu_u_sum-nu_u)
error2 = norm(nu_w_sum-nu_w)
