% make_data.m
%
% Script to run TFE with IIO and i*eta in polar (two layers)
%
% XT 5/18

clear all;
close all;
warning off;

% %% at 1-10^-16
% % 0.005
% a = 1-1e-16; 
% Eps = 0.005;
% setup_para;
% 
% tic;
% [I_u_n,I_w_n] = twolayer_IIO_TFE(zeta_n,psi_n,f,f_theta,p,k_u,k_w,sigma_u,sigma_w,a,b,c,N_theta,N,N_r,eta);
% [Un,Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,b,p,N_theta,N,N_r,sigma_u,eta);
% Qn_u = IIO_tfe_helmholtz_polar_exterior(Un,Dr_Un,Dp_Un,f,f_theta,a,b,N,sigma_u,eta);
% [Wn,Dr_Wn,Dp_Wn] = field_tfe_IIO_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,c,p,N_theta,N,N_r,sigma_w,eta);
% Sn_w = IIO_tfe_helmholtz_polar_interior(Wn,Dr_Wn,Dp_Wn,f,f_theta,a,c,N,sigma_w,eta);
% t_tfe = toc;
% 
% fprintf('  t_tfe = %g\n',t_tfe);
% 
% filename = sprintf('IIO_Eps_%g_Nr%g_eta%g_sing16.mat',Eps,N_r,eta);
% save(filename,'t_tfe','Eps','N','N_theta','N_r','lambda','k_u','k_w','a','b','c',...
%     'sigma_u','sigma_w','eta','I_u_n','I_w_n','Un','Qn_u','Wn','Sn_w',...
%     'xi_u','xi_w','Q_u','S_w');
% 
% % 0.01
% Eps = 0.01;
% setup_para;
% 
% tic;
% [I_u_n,I_w_n] = twolayer_IIO_TFE(zeta_n,psi_n,f,f_theta,p,k_u,k_w,sigma_u,sigma_w,a,b,c,N_theta,N,N_r,eta);
% [Un,Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,b,p,N_theta,N,N_r,sigma_u,eta);
% Qn_u = IIO_tfe_helmholtz_polar_exterior(Un,Dr_Un,Dp_Un,f,f_theta,a,b,N,sigma_u,eta);
% [Wn,Dr_Wn,Dp_Wn] = field_tfe_IIO_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,c,p,N_theta,N,N_r,sigma_w,eta);
% Sn_w = IIO_tfe_helmholtz_polar_interior(Wn,Dr_Wn,Dp_Wn,f,f_theta,a,c,N,sigma_w,eta);
% t_tfe = toc;
% 
% fprintf('  t_tfe = %g\n',t_tfe);
% 
% filename = sprintf('IIO_Eps_%g_Nr%g_eta%g_sing16.mat',Eps,N_r,eta);
% save(filename,'t_tfe','Eps','N','N_theta','N_r','lambda','k_u','k_w','a','b','c',...
%     'sigma_u','sigma_w','eta','I_u_n','I_w_n','Un','Qn_u','Wn','Sn_w',...
%     'xi_u','xi_w','Q_u','S_w');
% 
% % 0.05 
% Eps = 0.05;
% setup_para;
% 
% tic;
% [I_u_n,I_w_n] = twolayer_IIO_TFE(zeta_n,psi_n,f,f_theta,p,k_u,k_w,sigma_u,sigma_w,a,b,c,N_theta,N,N_r,eta);
% [Un,Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,b,p,N_theta,N,N_r,sigma_u,eta);
% Qn_u = IIO_tfe_helmholtz_polar_exterior(Un,Dr_Un,Dp_Un,f,f_theta,a,b,N,sigma_u,eta);
% [Wn,Dr_Wn,Dp_Wn] = field_tfe_IIO_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,c,p,N_theta,N,N_r,sigma_w,eta);
% Sn_w = IIO_tfe_helmholtz_polar_interior(Wn,Dr_Wn,Dp_Wn,f,f_theta,a,c,N,sigma_w,eta);
% t_tfe = toc;
% 
% fprintf('  t_tfe = %g\n',t_tfe);
% 
% filename = sprintf('IIO_Eps_%g_Nr%g_eta%g_sing16.mat',Eps,N_r,eta);
% save(filename,'t_tfe','Eps','N','N_theta','N_r','lambda','k_u','k_w','a','b','c',...
%     'sigma_u','sigma_w','eta','I_u_n','I_w_n','Un','Qn_u','Wn','Sn_w',...
%     'xi_u','xi_w','Q_u','S_w');
% 
% % 0.1
% Eps = 0.1;
% setup_para;
% 
% tic;
% [I_u_n,I_w_n] = twolayer_IIO_TFE(zeta_n,psi_n,f,f_theta,p,k_u,k_w,sigma_u,sigma_w,a,b,c,N_theta,N,N_r,eta);
% [Un,Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,b,p,N_theta,N,N_r,sigma_u,eta);
% Qn_u = IIO_tfe_helmholtz_polar_exterior(Un,Dr_Un,Dp_Un,f,f_theta,a,b,N,sigma_u,eta);
% [Wn,Dr_Wn,Dp_Wn] = field_tfe_IIO_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,c,p,N_theta,N,N_r,sigma_w,eta);
% Sn_w = IIO_tfe_helmholtz_polar_interior(Wn,Dr_Wn,Dp_Wn,f,f_theta,a,c,N,sigma_w,eta);
% t_tfe = toc;
% 
% fprintf('  t_tfe = %g\n',t_tfe);
% 
% filename = sprintf('IIO_Eps_%g_Nr%g_eta%g_sing16.mat',Eps,N_r,eta);
% save(filename,'t_tfe','Eps','N','N_theta','N_r','lambda','k_u','k_w','a','b','c',...
%     'sigma_u','sigma_w','eta','I_u_n','I_w_n','Un','Qn_u','Wn','Sn_w',...
%     'xi_u','xi_w','Q_u','S_w');

%% at 1-10^-12
% 0.005
% a = 1-1e-12  ; 
% Eps = 0.005;
% setup_para;
% 
% tic;
% [I_u_n,I_w_n] = twolayer_IIO_TFE(zeta_n,psi_n,f,f_theta,p,k_u,k_w,sigma_u,sigma_w,a,b,c,N_theta,N,N_r,eta);
% [Un,Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,b,p,N_theta,N,N_r,sigma_u,eta);
% Qn_u = IIO_tfe_helmholtz_polar_exterior(Un,Dr_Un,Dp_Un,f,f_theta,a,b,N,sigma_u,eta);
% [Wn,Dr_Wn,Dp_Wn] = field_tfe_IIO_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,c,p,N_theta,N,N_r,sigma_w,eta);
% Sn_w = IIO_tfe_helmholtz_polar_interior(Wn,Dr_Wn,Dp_Wn,f,f_theta,a,c,N,sigma_w,eta);
% t_tfe = toc;
% 
% fprintf('  t_tfe = %g\n',t_tfe);
% 
% filename = sprintf('IIO_Eps_%g_Nr%g_eta%g_sing12.mat',Eps,N_r,eta);
% save(filename,'t_tfe','Eps','N','N_theta','N_r','lambda','k_u','k_w','a','b','c',...
%     'sigma_u','sigma_w','eta','I_u_n','I_w_n','Un','Qn_u','Wn','Sn_w',...
%     'xi_u','xi_w','Q_u','S_w');
% 
% % 0.01
% Eps = 0.01;
% setup_para;
% 
% tic;
% [I_u_n,I_w_n] = twolayer_IIO_TFE(zeta_n,psi_n,f,f_theta,p,k_u,k_w,sigma_u,sigma_w,a,b,c,N_theta,N,N_r,eta);
% [Un,Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,b,p,N_theta,N,N_r,sigma_u,eta);
% Qn_u = IIO_tfe_helmholtz_polar_exterior(Un,Dr_Un,Dp_Un,f,f_theta,a,b,N,sigma_u,eta);
% [Wn,Dr_Wn,Dp_Wn] = field_tfe_IIO_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,c,p,N_theta,N,N_r,sigma_w,eta);
% Sn_w = IIO_tfe_helmholtz_polar_interior(Wn,Dr_Wn,Dp_Wn,f,f_theta,a,c,N,sigma_w,eta);
% t_tfe = toc;
% 
% fprintf('  t_tfe = %g\n',t_tfe);
% 
% filename = sprintf('IIO_Eps_%g_Nr%g_eta%g_sing12.mat',Eps,N_r,eta);
% save(filename,'t_tfe','Eps','N','N_theta','N_r','lambda','k_u','k_w','a','b','c',...
%     'sigma_u','sigma_w','eta','I_u_n','I_w_n','Un','Qn_u','Wn','Sn_w',...
%     'xi_u','xi_w','Q_u','S_w');
% 
% % 0.05 
% Eps = 0.05;
% setup_para;
% 
% tic;
% [I_u_n,I_w_n] = twolayer_IIO_TFE(zeta_n,psi_n,f,f_theta,p,k_u,k_w,sigma_u,sigma_w,a,b,c,N_theta,N,N_r,eta);
% [Un,Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,b,p,N_theta,N,N_r,sigma_u,eta);
% Qn_u = IIO_tfe_helmholtz_polar_exterior(Un,Dr_Un,Dp_Un,f,f_theta,a,b,N,sigma_u,eta);
% [Wn,Dr_Wn,Dp_Wn] = field_tfe_IIO_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,c,p,N_theta,N,N_r,sigma_w,eta);
% Sn_w = IIO_tfe_helmholtz_polar_interior(Wn,Dr_Wn,Dp_Wn,f,f_theta,a,c,N,sigma_w,eta);
% t_tfe = toc;
% 
% fprintf('  t_tfe = %g\n',t_tfe);
% 
% filename = sprintf('IIO_Eps_%g_Nr%g_eta%g_sing12.mat',Eps,N_r,eta);
% save(filename,'t_tfe','Eps','N','N_theta','N_r','lambda','k_u','k_w','a','b','c',...
%     'sigma_u','sigma_w','eta','I_u_n','I_w_n','Un','Qn_u','Wn','Sn_w',...
%     'xi_u','xi_w','Q_u','S_w');
% 
% % 0.1
% Eps = 0.1;
% setup_para;
% 
% tic;
% [I_u_n,I_w_n] = twolayer_IIO_TFE(zeta_n,psi_n,f,f_theta,p,k_u,k_w,sigma_u,sigma_w,a,b,c,N_theta,N,N_r,eta);
% [Un,Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,b,p,N_theta,N,N_r,sigma_u,eta);
% Qn_u = IIO_tfe_helmholtz_polar_exterior(Un,Dr_Un,Dp_Un,f,f_theta,a,b,N,sigma_u,eta);
% [Wn,Dr_Wn,Dp_Wn] = field_tfe_IIO_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,c,p,N_theta,N,N_r,sigma_w,eta);
% Sn_w = IIO_tfe_helmholtz_polar_interior(Wn,Dr_Wn,Dp_Wn,f,f_theta,a,c,N,sigma_w,eta);
% t_tfe = toc;
% 
% fprintf('  t_tfe = %g\n',t_tfe);
% 
% filename = sprintf('IIO_Eps_%g_Nr%g_eta%g_sing12.mat',Eps,N_r,eta);
% save(filename,'t_tfe','Eps','N','N_theta','N_r','lambda','k_u','k_w','a','b','c',...
%     'sigma_u','sigma_w','eta','I_u_n','I_w_n','Un','Qn_u','Wn','Sn_w',...
%     'xi_u','xi_w','Q_u','S_w');
% 
% 

%% at 0.5
% 0.005
a = 0.5;
% Eps = 0.005;
% setup_para;
% b = 0.8;
% c = 0.3;
% 
% tic;
% [I_u_n,I_w_n] = twolayer_IIO_TFE(zeta_n,psi_n,f,f_theta,p,k_u,k_w,sigma_u,sigma_w,a,b,c,N_theta,N,N_r,eta);
% [Un,Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,b,p,N_theta,N,N_r,sigma_u,eta);
% Qn_u = IIO_tfe_helmholtz_polar_exterior(Un,Dr_Un,Dp_Un,f,f_theta,a,b,N,sigma_u,eta);
% [Wn,Dr_Wn,Dp_Wn] = field_tfe_IIO_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,c,p,N_theta,N,N_r,sigma_w,eta);
% Sn_w = IIO_tfe_helmholtz_polar_interior(Wn,Dr_Wn,Dp_Wn,f,f_theta,a,c,N,sigma_w,eta);
% t_tfe = toc;
% 
% fprintf('  t_tfe = %g\n',t_tfe);
% 
% filename = sprintf('IIO_Eps_%g_Nr%g_eta%g.mat',Eps,N_r,eta);
% save(filename,'t_tfe','Eps','N','N_theta','N_r','lambda','k_u','k_w','a','b','c',...
%     'sigma_u','sigma_w','eta','I_u_n','I_w_n','Un','Qn_u','Wn','Sn_w',...
%     'xi_u','xi_w','Q_u','S_w');

% 0.01
Eps = 0.01;
setup_para;
b = 0.8;
c = 0.3;

tic;
[I_u_n,I_w_n] = twolayer_IIO_TFE(zeta_n,psi_n,f,f_theta,p,k_u,k_w,sigma_u,sigma_w,a,b,c,N_theta,N,N_r,eta);
[Un,Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,b,p,N_theta,N,N_r,sigma_u,eta);
Qn_u = IIO_tfe_helmholtz_polar_exterior(Un,Dr_Un,Dp_Un,f,f_theta,a,b,N,sigma_u,eta);
[Wn,Dr_Wn,Dp_Wn] = field_tfe_IIO_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,c,p,N_theta,N,N_r,sigma_w,eta);
Sn_w = IIO_tfe_helmholtz_polar_interior(Wn,Dr_Wn,Dp_Wn,f,f_theta,a,c,N,sigma_w,eta);
t_tfe = toc;

fprintf('  t_tfe = %g\n',t_tfe);

filename = sprintf('IIO_Eps_%g_Nr%g_eta%g.mat',Eps,N_r,eta);
save(filename,'t_tfe','Eps','N','N_theta','N_r','lambda','k_u','k_w','a','b','c',...
    'sigma_u','sigma_w','eta','I_u_n','I_w_n','Un','Qn_u','Wn','Sn_w',...
    'xi_u','xi_w','Q_u','S_w');

% 0.05 
Eps = 0.05;
setup_para;
b = 0.8;
c = 0.3;

tic;
[I_u_n,I_w_n] = twolayer_IIO_TFE(zeta_n,psi_n,f,f_theta,p,k_u,k_w,sigma_u,sigma_w,a,b,c,N_theta,N,N_r,eta);
[Un,Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,b,p,N_theta,N,N_r,sigma_u,eta);
Qn_u = IIO_tfe_helmholtz_polar_exterior(Un,Dr_Un,Dp_Un,f,f_theta,a,b,N,sigma_u,eta);
[Wn,Dr_Wn,Dp_Wn] = field_tfe_IIO_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,c,p,N_theta,N,N_r,sigma_w,eta);
Sn_w = IIO_tfe_helmholtz_polar_interior(Wn,Dr_Wn,Dp_Wn,f,f_theta,a,c,N,sigma_w,eta);
t_tfe = toc;

fprintf('  t_tfe = %g\n',t_tfe);

filename = sprintf('IIO_Eps_%g_Nr%g_eta%g.mat',Eps,N_r,eta);
save(filename,'t_tfe','Eps','N','N_theta','N_r','lambda','k_u','k_w','a','b','c',...
    'sigma_u','sigma_w','eta','I_u_n','I_w_n','Un','Qn_u','Wn','Sn_w',...
    'xi_u','xi_w','Q_u','S_w');

% 0.1
Eps = 0.1;
setup_para;
b = 0.8;
c = 0.3;

tic;
[I_u_n,I_w_n] = twolayer_IIO_TFE(zeta_n,psi_n,f,f_theta,p,k_u,k_w,sigma_u,sigma_w,a,b,c,N_theta,N,N_r,eta);
[Un,Dr_Un,Dp_Un] = field_tfe_IIO_helmholtz_polar_exterior(I_u_n,f,f_theta,k_u,a,b,p,N_theta,N,N_r,sigma_u,eta);
Qn_u = IIO_tfe_helmholtz_polar_exterior(Un,Dr_Un,Dp_Un,f,f_theta,a,b,N,sigma_u,eta);
[Wn,Dr_Wn,Dp_Wn] = field_tfe_IIO_helmholtz_polar_interior(I_w_n,f,f_theta,k_w,a,c,p,N_theta,N,N_r,sigma_w,eta);
Sn_w = IIO_tfe_helmholtz_polar_interior(Wn,Dr_Wn,Dp_Wn,f,f_theta,a,c,N,sigma_w,eta);
t_tfe = toc;

fprintf('  t_tfe = %g\n',t_tfe);

filename = sprintf('IIO_Eps_%g_Nr%g_eta%g.mat',Eps,N_r,eta);
save(filename,'t_tfe','Eps','N','N_theta','N_r','lambda','k_u','k_w','a','b','c',...
    'sigma_u','sigma_w','eta','I_u_n','I_w_n','Un','Qn_u','Wn','Sn_w',...
    'xi_u','xi_w','Q_u','S_w');



