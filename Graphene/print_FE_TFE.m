%% compare FE and TFE

clear all
close all

% load('IIO_fe_TM_Eps_40.mat');load('IIO_tfe_TM_Eps_40.mat'); %0.2
% load('IIO_fe_TM_Eps_20.mat');load('IIO_tfe_TM_Eps_20.mat'); %0.2
% load('IIO_fe_TM_Eps_1.mat');load('IIO_tfe_TM_Eps_1.mat'); %0.01
% load('IIO_fe_TM_cos4_Eps_40_Nr48.mat');load('IIO_tfe_TM_cos4_Eps_40_Nr48.mat'); %0.4
% load('IIO_fe_TM_cos4_Eps_20_Nr48.mat');load('IIO_tfe_TM_cos4_Eps_20_Nr48.mat'); 
% load('IIO_fe_TM_cos2_Eps_40_Nr48.mat');load('IIO_tfe_TM_cos2_Eps_40_Nr48.mat'); 

name = 'cos4'; 
% load('nIIO_fe_TM_cos4_Eps_1_Nr32_Nth128.mat');load('nIIO_tfe_TM_cos4_Eps_1_Nr32_Nth128.mat'); %0.01
load('nIIO_fe_TM_cos4_Eps_40_Nr32_Nth128.mat');load('nIIO_tfe_TM_cos4_Eps_40_Nr32_Nth128.mat'); %0.4
% load('nIIO_fe_TM_cos4_Eps_40_Nr48_Nth128.mat');load('nIIO_tfe_TM_cos4_Eps_40_Nr48_Nth128.mat'); %0.4
% load('nIIO_fe_TM_cos4_Eps_20_Nr32_Nth128.mat');load('nIIO_tfe_TM_cos4_Eps_20_Nr32_Nth128.mat'); %0.4

% name = 'cos2'; 
% load('nIIO_fe_TM_cos2_Eps_40_Nr32_Nth128.mat');load('nIIO_tfe_TM_cos4_Eps_40_Nr32_Nth128.mat'); %0.4
% load('nIIO_fe_TM_cos2_Eps_20_Nr32_Nth128.mat');load('nIIO_tfe_TM_cos4_Eps_20_Nr32_Nth128.mat'); %0.4



SavePlots = 1;
warning off
fprintf('test_FE_TFE_twolayer\n');
fprintf('-------------\n');
% fprintf('Lambda = %g  k_u = %g  k_w = %g\n',lambda,k_u,k_w);
% fprintf('Eps = %g  a = %g  b = %g  c = %g\n',Eps,a,b,c);
fprintf('N_theta = %d N = %d  N_r = %d\n',N_theta,N,N_r);
fprintf('\n');

[relerr,nplot] = compute_errors_2d_polar(S_w,S_w_n_fe,S_w_n_tfe,Eps,N,N_theta);
save_plots(SavePlots,nplot,relerr)


if(SavePlots==1)
  filename = sprintf('fig_err_fe_tfe_%s_Eps_%g_Nr%g_Nth%g',name,Eps*100/a,N_r,N_theta);
  saveas(1234,filename,'eps');
end

fprintf('\n');







