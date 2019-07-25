% test_IIO_FE_twolayer.m
%
% Script to test IIO in polar (two layers)
%
% XT 11/18

clear all;
close all;
warning off;




%% eta

%% at 1-10^-16
% 0.005
a = 1-1e-16; 
Eps = 0.005;
setup_para;

eta = 3.4;
Z_p = -1i*eta.*ones(N_theta,1);
Y_p = 1i*eta.*ones(N_theta,1);
I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_fe_Eps_%g_eta%g_sing16.mat',Eps,eta);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','eta','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');

% 0.01
Eps = 0.01;
setup_para;

eta = 3.4;
Z_p = -1i*eta.*ones(N_theta,1);
Y_p = 1i*eta.*ones(N_theta,1);
I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_fe_Eps_%g_eta%g_sing16.mat',Eps,eta);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','eta','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');

% 0.05 
Eps = 0.05;
setup_para;

eta = 3.4;
Z_p = -1i*eta.*ones(N_theta,1);
Y_p = 1i*eta.*ones(N_theta,1);
I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_fe_Eps_%g_eta%g_sing16.mat',Eps,eta);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','eta','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');

% 0.1
Eps = 0.1;
setup_para;

eta = 3.4;
Z_p = -1i*eta.*ones(N_theta,1);
Y_p = 1i*eta.*ones(N_theta,1);
I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_fe_Eps_%g_eta%g_sing16.mat',Eps,eta);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','eta','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');


% at 1-10^-12
% 0.005
a = 1-1e-12; 
Eps = 0.005;
setup_para;

eta = 3.4;
Z_p = -1i*eta.*ones(N_theta,1);
Y_p = 1i*eta.*ones(N_theta,1);
I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_fe_Eps_%g_eta%g_sing12.mat',Eps,eta);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','eta','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');

% 0.01
Eps = 0.01;
setup_para;

eta = 3.4;
Z_p = -1i*eta.*ones(N_theta,1);
Y_p = 1i*eta.*ones(N_theta,1);
I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_fe_Eps_%g_eta%g_sing12.mat',Eps,eta);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','eta','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');

% 0.05 
Eps = 0.05;
setup_para;

eta = 3.4;
Z_p = -1i*eta.*ones(N_theta,1);
Y_p = 1i*eta.*ones(N_theta,1);
I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_fe_Eps_%g_eta%g_sing12.mat',Eps,eta);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','eta','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');

% 0.1
Eps = 0.1;
setup_para;

eta = 3.4;
Z_p = -1i*eta.*ones(N_theta,1);
Y_p = 1i*eta.*ones(N_theta,1);
I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_fe_Eps_%g_eta%g_sing12.mat',Eps,eta);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','eta','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');

% at 0.5
% 0.005
a = 0.5;
Eps = 0.005;
setup_para;

eta = 3.4;
Z_p = -1i*eta.*ones(N_theta,1);
Y_p = 1i*eta.*ones(N_theta,1);
I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_fe_Eps_%g_eta%g.mat',Eps,eta);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','eta','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');

% 0.01
Eps = 0.01;
setup_para;

eta = 3.4;
Z_p = -1i*eta.*ones(N_theta,1);
Y_p = 1i*eta.*ones(N_theta,1);
I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_fe_Eps_%g_eta%g.mat',Eps,eta);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','eta','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');

% 0.05 
Eps = 0.05;
setup_para;

eta = 3.4;
Z_p = -1i*eta.*ones(N_theta,1);
Y_p = 1i*eta.*ones(N_theta,1);
I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_fe_Eps_%g_eta%g.mat',Eps,eta);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','eta','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');

% 0.1
Eps = 0.1;
setup_para;

eta = 3.4;
Z_p = -1i*eta.*ones(N_theta,1);
Y_p = 1i*eta.*ones(N_theta,1);
I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_fe_Eps_%g_eta%g.mat',Eps,eta);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','eta','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');


%% S0=Q0=0

%% at 1-10^-16
% 0.005
a = 1-1e-16; 
Eps = 0.005;
setup_para;
Z_p = sigma_u * k_u * diff_besselh(p,1,k_u*a)./besselh(p,k_u*a);
Y_p = sigma_w * k_w * diff_besselj(p,1,k_w*a)./besselj(p,k_w*a);

I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_new_fe_Eps_%g_sing16.mat',Eps);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','Y_p','Z_p','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');

% 0.01
Eps = 0.01;
setup_para;
Z_p = sigma_u * k_u * diff_besselh(p,1,k_u*a)./besselh(p,k_u*a);
Y_p = sigma_w * k_w * diff_besselj(p,1,k_w*a)./besselj(p,k_w*a);

I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_new_fe_Eps_%g_sing16.mat',Eps);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','Y_p','Z_p','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');

% 0.05 
Eps = 0.05;
setup_para;
Z_p = sigma_u * k_u * diff_besselh(p,1,k_u*a)./besselh(p,k_u*a);
Y_p = sigma_w * k_w * diff_besselj(p,1,k_w*a)./besselj(p,k_w*a);

I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_new_fe_Eps_%g_sing16.mat',Eps);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','Y_p','Z_p','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');

% 0.1
Eps = 0.1;
setup_para;
Z_p = sigma_u * k_u * diff_besselh(p,1,k_u*a)./besselh(p,k_u*a);
Y_p = sigma_w * k_w * diff_besselj(p,1,k_w*a)./besselj(p,k_w*a);

I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_new_fe_Eps_%g_sing16.mat',Eps);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','Y_p','Z_p','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');

% at 1-10^-12
% 0.005
a = 1-1e-12; 
Eps = 0.005;
setup_para;
Z_p = sigma_u * k_u * diff_besselh(p,1,k_u*a)./besselh(p,k_u*a);
Y_p = sigma_w * k_w * diff_besselj(p,1,k_w*a)./besselj(p,k_w*a);

I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_new_fe_Eps_%g_sing12.mat',Eps);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','Y_p','Z_p','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');

% 0.01
Eps = 0.01;
setup_para;
Z_p = sigma_u * k_u * diff_besselh(p,1,k_u*a)./besselh(p,k_u*a);
Y_p = sigma_w * k_w * diff_besselj(p,1,k_w*a)./besselj(p,k_w*a);

I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_new_fe_Eps_%g_sing12.mat',Eps);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','Y_p','Z_p','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');

% 0.05 
Eps = 0.05;
setup_para;
Z_p = sigma_u * k_u * diff_besselh(p,1,k_u*a)./besselh(p,k_u*a);
Y_p = sigma_w * k_w * diff_besselj(p,1,k_w*a)./besselj(p,k_w*a);

I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_new_fe_Eps_%g_sing12.mat',Eps);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','Y_p','Z_p','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');

% 0.1
Eps = 0.1;
setup_para;
Z_p = sigma_u * k_u * diff_besselh(p,1,k_u*a)./besselh(p,k_u*a);
Y_p = sigma_w * k_w * diff_besselj(p,1,k_w*a)./besselj(p,k_w*a);

I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_new_fe_Eps_%g_sing12.mat',Eps);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','Y_p','Z_p','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');



% at 0.5
% 0.005
a = 0.5;
Eps = 0.005;
setup_para;
Z_p = sigma_u * k_u * diff_besselh(p,1,k_u*a)./besselh(p,k_u*a);
Y_p = sigma_w * k_w * diff_besselj(p,1,k_w*a)./besselj(p,k_w*a);

I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_new_fe_Eps_%g.mat',Eps);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','Y_p','Z_p','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');

% 0.01
Eps = 0.01;
setup_para;
Z_p = sigma_u * k_u * diff_besselh(p,1,k_u*a)./besselh(p,k_u*a);
Y_p = sigma_w * k_w * diff_besselj(p,1,k_w*a)./besselj(p,k_w*a);

I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_new_fe_Eps_%g.mat',Eps);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','Y_p','Z_p','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');

% 0.05 
Eps = 0.05;
setup_para;
Z_p = sigma_u * k_u * diff_besselh(p,1,k_u*a)./besselh(p,k_u*a);
Y_p = sigma_w * k_w * diff_besselj(p,1,k_w*a)./besselj(p,k_w*a);

I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_new_fe_Eps_%g.mat',Eps);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','Y_p','Z_p','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');

% 0.1
Eps = 0.1;
setup_para;
Z_p = sigma_u * k_u * diff_besselh(p,1,k_u*a)./besselh(p,k_u*a);
Y_p = sigma_w * k_w * diff_besselj(p,1,k_w*a)./besselj(p,k_w*a);

I_u = sigma_u*nu_u+ifft(Y_p.*fft(xi_u));
I_w = sigma_w*nu_w-ifft(Z_p.*fft(xi_w));
Q_u = sigma_u*nu_u+ifft(Z_p.*fft(xi_u));
S_w = sigma_w*nu_w-ifft(Y_p.*fft(xi_w));

tic;
[I_u_n,I_w_n] = twolayer_IIO_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,...
    p,k_u,k_w,sigma_u,sigma_w,a,N_theta,N,Y_p,Z_p);
anp = field_fe_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Y_p);
Q_u_n = IIO_fe_helmholtz_polar_exterior(anp,f,f_theta,k_u,a,p,N_theta,N,sigma_u,Z_p);
dnp = field_fe_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Z_p);
S_w_n = IIO_fe_helmholtz_polar_interior(dnp,f,f_theta,k_w,a,p,N_theta,N,sigma_w,Y_p);
t_fe = toc;

fprintf('  t_fe = %g\n',t_fe);

filename = sprintf('IIO_new_fe_Eps_%g.mat',Eps);
save(filename,'t_fe','Eps','N','N_theta','lambda','k_u','k_w','a',...
    'sigma_u','sigma_w','Y_p','Z_p','I_u_n','I_w_n','Q_u_n','S_w_n',...
    'I_u','I_w','Q_u','S_w');
