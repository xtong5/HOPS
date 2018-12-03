
% test_helmholtz_polar_twolayer.m
%
% Script to test Helmholtz DNO solvers in polar (two layers)
%
% XT 12/16
% XT 11/18 save data

clear all;
close all;
warning off;
SavePlots = 0;
SaveData = 1;

RunNumber = 1;
Mode = 2; %check 


L = 2*pi;
lambda = 0.45;
n_u = 1;
% n_w = 2.5;
k_0 = L/lambda;
k_u = n_u*k_0; 
% k_w = n_w*k_0;
k_w = 5.13562230184068;

if(RunNumber==1)
  % Small Deformation
  Eps = 0.1;
  N_theta = 64;
%   a = 0.5;
%   a = 1-1e-16;
  a = 1-1e-12;
%   b = 10*a;
  N = 16;
elseif(RunNumber==2)
  % Big Deformation (inside disk)
  Eps = 0.3;
  N_theta = 64;
  a = 2.0;
  N = 16;
elseif(RunNumber==3)
  % Big Deformation (outside disk)
  Eps = 0.75;
  N_theta = 64;
  a = 2.0;
  N = 16;
end

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

Ar_u = 2; pp = 2; % take a special wavenumber
Ar_w = 1; %r = 2; % take a special wavenumber ???same

xi_u = Ar_u*besselh(pp,k_u.*A).*exp(1i*pp.*theta);
nu_u = Ar_u*((-k_u.*A.*(diff_bessel(2,pp,1,k_u.*A))+...
    1i*pp*Eps.*f_theta.*besselh(pp,k_u.*A)./A ).*exp(1i*pp.*theta)); 
xi_w = Ar_w*besselj(pp,k_w.*A).*exp(1i*pp.*theta);
nu_w = Ar_w*((k_w.*A.*(diff_bessel(1,pp,1,k_w.*A))-...
    1i*pp*Eps.*f_theta.*besselj(pp,k_w.*A)./A ).*exp(1i*pp.*theta)); 
xi_u_n = zeros(N_theta,N+1); nu_u_n = zeros(N_theta,N+1);
xi_w_n = zeros(N_theta,N+1); nu_w_n = zeros(N_theta,N+1);
f_n = ones(N_theta,1); f_nmo = ones(N_theta,1);f_nmt = ones(N_theta,1);
f_n = f.*f_n;
xi_u_n(:,0+1) = Ar_u*besselh(pp,k_u*a).*exp(1i*pp.*theta);
xi_u_n(:,1+1) = Ar_u*k_u^1*diff_bessel(2,pp,1,k_u*a).*f_n.*exp(1i*pp.*theta);
nu_u_n(:,0+1) = -Ar_u*k_u*a*diff_bessel(2,pp,1,k_u*a).*exp(1i*pp.*theta);
nu_u_n(:,1+1) = -f/a.*nu_u_n(:,1)...
      -Ar_u*a*k_u^(1+1).*diff_bessel(2,pp,1+1,k_u*a).*f_n.*exp(1i*pp.*theta)...
      -Ar_u*(2*f).*k_u^1.*diff_bessel(2,pp,1,k_u*a).*f_nmo.*exp(1i*pp.*theta)...
      +Ar_u*(f_theta/a).*(1i*pp).*besselh(pp,k_u*a).*f_nmo.*exp(1i*pp.*theta);
xi_w_n(:,0+1) = Ar_w*besselj(pp,k_w*a).*exp(1i*pp.*theta);
xi_w_n(:,1+1) = Ar_w*k_w^1*diff_bessel(1,pp,1,k_w*a).*f_n.*exp(1i*pp.*theta);
nu_w_n(:,0+1) = Ar_w*k_w*a*diff_bessel(1,pp,1,k_w*a).*exp(1i*pp.*theta);
nu_w_n(:,1+1) = -f/a.*nu_w_n(:,1)...
      +Ar_w*a*k_w^(1+1).*diff_bessel(1,pp,1+1,k_w*a).*f_n.*exp(1i*pp.*theta)...
      +Ar_w*(2*f).*k_w^1.*diff_bessel(1,pp,1,k_w*a).*f_nmo.*exp(1i*pp.*theta)...
      -Ar_w*(f_theta/a).*(1i*pp).*besselj(pp,k_w*a).*f_nmo.*exp(1i*pp.*theta);

for n=2:N
  f_n = f.*f_n/n;
  f_nmo = f.*f_nmo/(n-1);
  if(n>2)
    f_nmt = f.*f_nmt/(n-2);
  end
  xi_u_n(:,n+1) = Ar_u*k_u^n*diff_bessel(2,pp,n,k_u*a).*f_n.*exp(1i*pp.*theta);
  nu_u_n(:,n+1) = -f/a.*nu_u_n(:,n-1+1)...
      -Ar_u*a*k_u^(n+1).*diff_bessel(2,pp,n+1,k_u*a).*f_n.*exp(1i*pp.*theta)...
      -Ar_u*(2*f).*k_u^n.*diff_bessel(2,pp,n,k_u*a).*f_nmo.*exp(1i*pp.*theta)...
      -Ar_u*(f.^2/a)*k_u^(n-1).*diff_bessel(2,pp,n-1,k_u*a).*f_nmt.*exp(1i*pp.*theta)...
      +Ar_u*(f_theta/a)*k_u^(n-1).*(1i*pp).*diff_bessel(2,pp,n-1,k_u*a)...
      .*f_nmo.*exp(1i*pp.*theta);
  xi_w_n(:,n+1) = Ar_w*k_w^n*diff_bessel(1,pp,n,k_w*a).*f_n.*exp(1i*pp.*theta);
  nu_w_n(:,n+1) = -f/a.*nu_w_n(:,n-1+1)...
      +Ar_w*a*k_w^(n+1).*diff_bessel(1,pp,n+1,k_w*a).*f_n.*exp(1i*pp.*theta)...
      +Ar_w*(2*f).*k_w^n.*diff_bessel(1,pp,n,k_w*a).*f_nmo.*exp(1i*pp.*theta)...
      +Ar_w*(f.^2/a)*k_w^(n-1).*diff_bessel(1,pp,n-1,k_w*a).*f_nmt.*exp(1i*pp.*theta)...
      -Ar_w*(f_theta/a)*k_w^(n-1).*(1i*pp).*diff_bessel(1,pp,n-1,k_w*a)...
      .*f_nmo.*exp(1i*pp.*theta);
end


if(Mode==1)
  tau2 = 1;
else
  tau2 = k_u^2/k_w^2;
end
zeta_n = xi_u_n - xi_w_n; % nu_u points downwards!
psi_n = -nu_u_n - tau2*nu_w_n;



% Two-layer scattering by DNO


fprintf('\n\nTwo-layer scattering by DNO\n\n');

tic;
U_n = twolayer_dno_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,tau2,...
    p,k_u,k_w,a,N_theta,N);
apn_fe = field_fe_helmholtz_polar_exterior(U_n,f,k_u,a,p,N_theta,N);
Gn_fe_u = dno_fe_helmholtz_polar_exterior(apn_fe,f,f_theta,k_u,a,p,N_theta,N);
W_n = U_n - zeta_n;
dpn_fe = field_fe_helmholtz_polar_interior(W_n,f,k_w,a,p,N_theta,N);
Gn_fe_w = dno_fe_helmholtz_polar_interior(dpn_fe,f,f_theta,k_w,a,p,N_theta,N);
t_fe = toc;

if SaveData==1
% filename = sprintf('DNO_fe_Eps_%g.mat',Eps);
filename = sprintf('DNO_fe_Eps_%g_sing12.mat',Eps);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'U_n','W_n','Gn_fe_u','Gn_fe_w','xi_u','xi_w','nu_u','nu_w');
end

% fprintf('Press key to compute exterior layer errors...\n');
% pause;

fprintf('  t_fe = %g\n',t_fe);
% fprintf('\nEXTERIOR LAYER\n\n');
[relerrU,nplotU] = compute_errors_2d_polar(xi_u,U_n,Eps,N,N_theta);
[relerrDNOU,nplotDNOU] = compute_errors_2d_polar(nu_u,Gn_fe_u,Eps,N,N_theta);
% make_plots_polar(SavePlots,nplotU,relerrU);
% make_plots_polar(SavePlots,nplotDNOU,relerrDNOU);
% fprintf('\n');

% fprintf('Press key to compute interior layer errors...\n');
% % pause;
% 
% fprintf('\nINTERIOR LAYER\n\n');
% % [relerrW,nplotW] = compute_errors_2d_polar(xi_w,W_n,Eps,N,N_theta);
% % [relerrDNOW,nplotDNOW] = compute_errors_2d_polar(nu_w,Gn_fe_w,Eps,N,N_theta);
% % make_plots_polar(SavePlots,nplotW,relerrW);
% % make_plots_polar(SavePlots,nplotDNOW,relerrDNOW);
% fprintf('\n');

% fprintf('Press key to compute the far field behavior...\n');
% pause;

% hold on;
% bb = [1 5 10 20 40 80];
% for jj=1:length(bb)
%   b = bb(jj);
%   B_far = zeros(N_theta,1);
%   coeff=zeros(N_theta,N+1);
%   for n=1:N+1
%     coeff(:,n)=ifft(apn_fe(:,n).*besselh(p,k_u.*b)./besselh(p,k_u.*a)); 
%   end 
% 
%   for j=1:N_theta
%     B_far(j) = taylorsum(coeff(j,:).',Eps,N);
%   end
%   if(jj>1)
%     u_tilde_prev = u_tilde;
%   end
%   u_tilde= sqrt(b)*exp(-1i*k_u*b).*B_far;
%   if(jj>1)
%     fprintf('  |u_tilde-u_tilde_prev| = %g\n',norm(u_tilde-u_tilde_prev,inf));
%   end
% %   
%   
%   if(jj==1)
%     plot(theta,real(u_tilde),'b-o');
%   elseif(jj==2)
%     plot(theta,real(u_tilde),'r-o');
%   elseif(jj==3)
%     plot(theta,real(u_tilde),'g-o');
%   elseif(jj==4)
%     plot(theta,real(u_tilde),'c-o');
%   elseif(jj==5)
%     plot(theta,real(u_tilde),'y-o');
%   else
%     plot(theta,real(u_tilde),'k-o');
%   end
%   hold on;
% end
% 
% C_b = sqrt(b)*exp(-1i*k_u*b);
% B_far_true = Ar_u*besselh(r,k_u.*b).*exp(1i*r.*theta);
% B_far = zeros(N_theta,N+1);
% B_far1 = zeros(N_theta,1);
% for n=1:N+1
%   B_far(:,n)=ifft(apn_fe(:,n).*besselh(p,k_u.*b)./besselh(p,k_u.*a)); 
% end 
% 
% for j=1:N_theta
%   B_far1(j) = taylorsum(B_far(j,:),Eps,N);
% end

% [relerrB,nplotB] = compute_errors_2d_polar(C_b*B_far_true,C_b*B_far,Eps,N,N_theta);
% make_plots_polar(SavePlots,nplotB,relerrB);
% plot(1:N_theta,real(B_far_true),'b-o',1:N_theta,real(B_far1),'r-*')

% bb = [1 5 10 20 40 80];
% for jj=1:length(bb)
%   bbb = bb(jj);
%   BB_far = zeros(N_theta,1);
%   Bcoeff=zeros(N_theta,N+1);
%   for n=1:N+1
%     Bcoeff(:,n)=ifft(apn_fe(:,n).*besselh(p,k_u.*bbb)./besselh(p,k_u.*a)); 
%   end 
%   for j=1:N_theta
%     B_far(j) = taylorsum(Bcoeff(j,:).',Eps,N);
%   end
%   if(jj>1)
%     u_tilde_prev = u_tilde;
%   end
%   u_tilde= sqrt(bbb)*exp(-1i*k_u*bbb).*BB_far;
%   if(jj>1)
%   BB_error = log(abs(u_tilde-u_tilde_prev));
%   end
% end