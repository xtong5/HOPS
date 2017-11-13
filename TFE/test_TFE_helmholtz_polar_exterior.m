clear all


%% Default
SavePlots = 0;
RunNumber = 1;
L = 2*pi;
k=2.1;

if(RunNumber==1)
  % Small Deformation
  Eps = 0.02;
  N_theta = 64;
  N = 16;
  N_r = 16;
  a = 1;
  b = 1.6;
elseif(RunNumber==2)
  % Big Deformation (inside disk)
  Eps = 0.3;
  N_theta = 64;
  N = 16;
  N_r = 16;
  a = 1;
  b = 1.6;
elseif(RunNumber==3)
  % Big Deformation (outside disk)
  Eps = 0.75;
  N_theta = 64;
  N = 16;
  N_r = 16;
  a = 1;
  b = 1.6;
end

p = [0:N_theta/2-1,-N_theta/2:-1]';
theta = (L/N_theta)*[0:N_theta-1]';

fprintf('test_TFE_helmholtz_polar_exterior\n');
fprintf('-------------\n');
fprintf('k = %g a = %g b = %g\n',k,a,b);
fprintf('Eps = %g \n',Eps);
fprintf('N_theta = %d N_r = %d N = %d\n',N_theta,N_r,N);
fprintf('\n');


%% data
f = exp(cos(theta));
f_theta = ifft( (1i*p).*fft(f) );

%% Exact Solution
Ar = 2; pp = 2; % compute a special wavenumber
xi = Ar*besselh(pp,k.*(a+Eps.*f)).*exp(1i*pp.*theta);
% u_exact = Ar * kron(exp(1i*pp*theta),besselh(pp,k*r));
u_exact = Ar * exp(1i*pp*theta) .* besselh(pp,k.*(a+Eps.*f));
% u_exact = Ar * exp(1i*pp*theta) .* besselh(pp,k.*(b));
% xi_n_p = xi_n_exterior(N_theta,N,f,Ar,pp,a,k,theta);
f_n = ones(N_theta,1); 
xi_n(:,0+1) = Ar*besselh(pp,k*a).*exp(1i*pp.*theta);
for n=1:N
  f_n = f.*f_n/n;
  xi_n(:,n+1) = Ar*k^n*diff_besselh(pp,n,k*a).*f_n.*exp(1i*pp.*theta);
end
nu =Ar*(-k.*(a+Eps.*f).*(diff_besselh(pp,1,k.*(a+Eps.*f)))+1i*pp*Eps.*...
    f_theta.*besselh(pp,k.*(a+Eps.*f))./(a+Eps.*f) ).*exp(1i*pp.*theta); %DNO


[Un,Dr_Un,Dp_Un] = field_tfe_helmholtz_polar_exterior(xi_n,f,f_theta,k,a,b,p,N_theta,N,N_r);
% [relerr,nplot] = compute_errors_2d_polar(u_exact,Un,Eps,N,N_theta);

Gn_tfe = dno_tfe_helmholtz_polar_exterior(Dr_Un,Dp_Un,f,f_theta,k,a,b,p,N_theta,N,N_r);
[relerr,nplot] = compute_errors_2d_polar(nu,Gn_tfe,Eps,N,N_theta);        


