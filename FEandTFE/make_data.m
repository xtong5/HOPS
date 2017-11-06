clear all;
close all;


Mode = 2; %check 
% Mode = 1;
L = 2*pi;
lambda = 0.45;
n_u = 1;
n_w = 2.5;
k_0 = L/lambda;
k_u = n_u*k_0; 
k_w = n_w*k_0;

%% first
Eps = 0.002;
N_theta = 64;
a = 0.025;
b = 10*a;
c = 0.1*a;
N = 16;
N_r = 32;
theta = (L/N_theta)*[0:N_theta-1]';
p = [0:N_theta/2-1,-N_theta/2:-1]';
f = exp(cos(theta));
f_theta = ifft( (1i*p).*fft(f) );
DirichletNuemannData
if(Mode==1)
  tau2 = 1;
else
  tau2 = k_u^2/k_w^2;
end
zeta_n = xi_u_n - xi_w_n;
% nu_u points downwards!
psi_n = -nu_u_n - tau2*nu_w_n;
%psi_n = nu_u_n - tau2*nu_w_n;

%% Two-layer scattering by DNO
name='Runnumber1';
filename = sprintf('data_%s.mat',name);
FE_TFE_twolayer(Eps,lambda,p,k_u,k_w,N_theta,N,N_r,a,b,c,f,f_theta,zeta_n,psi_n,tau2,filename);

%% second run
Eps = 0.002;
N_theta = 64;
a = 0.025;
b = 10*a;
c = 0.5*a;
N = 16;
N_r = 32; 
theta = (L/N_theta)*[0:N_theta-1]';
p = [0:N_theta/2-1,-N_theta/2:-1]';
f = exp(cos(theta));
f_theta = ifft( (1i*p).*fft(f) );
DirichletNuemannData
if(Mode==1)
  tau2 = 1;
else
  tau2 = k_u^2/k_w^2;
end
zeta_n = xi_u_n - xi_w_n;
% nu_u points downwards!
psi_n = -nu_u_n - tau2*nu_w_n;
%psi_n = nu_u_n - tau2*nu_w_n;

% Two-layer scattering by DNO
name='Runnumber2';
filename = sprintf('data_%s.mat',name);
FE_TFE_twolayer(Eps,lambda,p,k_u,k_w,N_theta,N,N_r,a,b,c,f,f_theta,zeta_n,psi_n,tau2,filename);


%% third
Eps = 2;
N_theta = 64;
a = 0.025;
b = 10*a;
c = 0.1*a;
N = 16;
N_r = 32;  
theta = (L/N_theta)*[0:N_theta-1]';
p = [0:N_theta/2-1,-N_theta/2:-1]';
f = exp(cos(theta));
f_theta = ifft( (1i*p).*fft(f) );
DirichletNuemannData
if(Mode==1)
  tau2 = 1;
else
  tau2 = k_u^2/k_w^2;
end
zeta_n = xi_u_n - xi_w_n;
% nu_u points downwards!
psi_n = -nu_u_n - tau2*nu_w_n;
%psi_n = nu_u_n - tau2*nu_w_n;

%% Two-layer scattering by DNO
name='Runnumber3';
filename = sprintf('data_%s.mat',name);
FE_TFE_twolayer(Eps,lambda,p,k_u,k_w,N_theta,N,N_r,a,b,c,f,f_theta,zeta_n,psi_n,tau2,filename);

