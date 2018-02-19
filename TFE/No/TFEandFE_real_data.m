close all
SavePlots = 0;
warning off


load('TFE_cos4_eps_a5_Nr64_Nt128_WATERAg.mat') % paper
load('FE_cos4_eps_a5_Nt128_WATERAg1.mat')
% load('TFE_cos8_eps_a10_Nr64_N16_WATERAg.mat') % paper
% load('FE_cos8_eps_a10_N16_WATERAg.mat')
% load('TFE_cos8_eps_a10_Nr64_N24_WATERAg.mat') % paper
% load('FE_cos8_eps_a10_N24_WATERAg.mat')



fprintf('lambda = %g  k_u = %g  k_w = %g\n\n',lambda,k_u,k_w);
fprintf('Eps = %g  a = %g  b = %g  c = %g\n',Eps,a,b,c);
fprintf('N_theta = %d N = %d  N_r = %d\n',N_theta,N,N_r);
fprintf('\n');


Gn_tfe =zeros(N_theta,1);
Gn_fe =zeros(N_theta,1);

conv_err_fe = zeros(5,1); %diff between each order 
conv_err_tfe = zeros(5,1);

i=1;M=[0 2.^[1:4]];
% i=1;M=[0 2.^[1:4] 20 24];
for i=1:size(M,2)
    for j=1:N_theta
        m = floor(M(i)/2);
        Gn_tfe(j) = padesum(Gn_tfe_u(j,:).',Eps,m);  
        Gn_fe(j) = padesum(Gn_fe_u(j,:).',Eps,m);  
    end
    Gn_U_tfe_norm(i) = norm(Gn_tfe)/sqrt(N_theta);
    Gn_U_fe_norm(i) = norm(Gn_fe)/sqrt(N_theta);
end

for i=2:size(M,2)
    conv_err_fe(i) = abs(Gn_U_fe_norm(i)-Gn_U_fe_norm(i-1));
    conv_err_tfe(i) = abs(Gn_U_tfe_norm(i)-Gn_U_tfe_norm(i-1));
end
    

fprintf('n      FE(P)      Err      TFE(P)      Err      Diff \n');
fprintf('-----------------------------------------------------\n');
for i=1:size(M,2)
  fprintf('%d    %.16g    %g      %.16g      %g      %g\n',M(i),...
  Gn_U_fe_norm(i),conv_err_fe(i),Gn_U_tfe_norm(i),conv_err_tfe(i),...
      abs(Gn_U_fe_norm(i)-Gn_U_tfe_norm(i)));
end






