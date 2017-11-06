%% print the results 

load('data_Runnumber1')
% load('data_Runnumber2')
% load('data_Runnumber3')
fprintf('test_FE_TFE_twolayer\n');
fprintf('-------------\n');
fprintf('Lambda = %g  k_u = %g  k_w = %g\n',lambda,k_u,k_w);
fprintf('Eps = %g  a = %g  b = %g  c = %g\n',Eps,a,b,c);
fprintf('N_theta = %d N = %d  N_r = %d\n',N_theta,N,N_r);
fprintf('\n');

fprintf('Press key to compute exterior layer errors...\n');
pause;

fprintf('  t_fe = %g  t_tfe = %g\n',t_fe,t_tfe);

fprintf('\nEXTERIOR LAYER\n\n');
[relerrDNOU,nplotDNOU] = compute_errors_2d_polar(nu_u,Gn_fe_u,Gn_tfe_u,Eps,N,N_theta);
% make_plots_polar(SavePlots,nplotDNOU,relerrDNOU);
fprintf('\n');

fprintf('Press key to compute interior layer errors...\n');
pause;
fprintf('\nINTERIOR LAYER\n\n');
[relerrDNOW,nplotDNOW] = compute_errors_2d_polar(nu_w,Gn_fe_w,Gn_tfe_w,Eps,N,N_theta);
% make_plots_polar(SavePlots,nplotDNOW,relerrDNOW);





