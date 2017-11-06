function [U_n_fe,W_n_fe,U_n_tfe,W_n_tfe,Gn_fe_u,Gn_fe_w,Gn_tfe_u,Gn_tfe_w]...
    =FE_TFE_twolayer(Eps,lambda,p,k_u,k_w,N_theta,N,N_r,a,b,c,f,f_theta,zeta_n,psi_n,tau2,filename)

% Two-layer scattering by DNO

%% FE
tic;
U_n_fe = twolayer_dno_fe_helmholtz_polar(zeta_n,psi_n,f,f_theta,tau2,...
    p,k_u,k_w,a,N_theta,N);
apn_fe = field_fe_helmholtz_polar_exterior(U_n_fe,f,k_u,a,p,N_theta,N);
Gn_fe_u = dno_fe_helmholtz_polar_exterior(apn_fe,f,f_theta,k_u,a,p,N_theta,N);
W_n_fe = U_n_fe - zeta_n;
dpn_fe = field_fe_helmholtz_polar_interior(W_n_fe,f,k_w,a,p,N_theta,N);
Gn_fe_w = dno_fe_helmholtz_polar_interior(dpn_fe,f,f_theta,k_w,a,p,N_theta,N);
t_fe = toc;

%% TFE
tic;
U_n_tfe = twolayer_dno_tfe_helmholtz_polar(zeta_n,psi_n,f,f_theta,tau2,...
    p,k_u,k_w,a,b,c,N_theta,N,N_r);
[Un,Dr_Un,Dp_Un] = field_tfe_helmholtz_polar_exterior(U_n_tfe,f,f_theta,k_u,a,b,p,N_theta,N,N_r);
Gn_tfe_u = dno_tfe_helmholtz_polar_exterior(Dr_Un,Dp_Un,f,f_theta,k_u,a,b,p,N_theta,N,N_r);
W_n_tfe = U_n_tfe - zeta_n;
[Wn,Dr_Wn,Dp_Wn] = field_tfe_helmholtz_polar_interior(W_n_tfe,f,f_theta,k_w,a,c,p,N_theta,N,N_r);
Gn_tfe_w = dno_tfe_helmholtz_polar_interior(Dr_Wn,Dp_Wn,f,f_theta,k_w,a,c,p,N_theta,N,N_r);
t_tfe = toc;

%% Save data
save(filename,'Eps','N','N_theta','N_r','lambda','k_u','k_w','a','b','c','tau2',...
    't_fe','t_tfe','U_n_fe','W_n_fe','U_n_tfe','W_n_tfe','Gn_fe_u','Gn_fe_w','Gn_tfe_u','Gn_tfe_w');

end