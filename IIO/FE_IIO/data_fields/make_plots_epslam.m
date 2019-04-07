% plot
clear all
close all

SavePlots=0;
function_mode = 2; % choose functions

if function_mode == 1  
    name = 'expcos';
%     load('data_expcos_eps10_VACAg.mat');   
    load('data_expcos_eps20_VACAg.mat'); 
%     load('data_expcos_eps10_WATERAg.mat');   
%     load('data_expcos_eps20_WATERAg.mat');   
end
if function_mode == 2
    name = 'cos2';
%     load('data_cos2_eps10_VACAg.mat');   
    load('data_cos2_eps20_VACAg.mat');   
%     load('data_cos2_eps10_WATERAg.mat');   
%     load('data_cos2_eps20_WATERAg.mat');   
end
if function_mode == 4
    name = 'cos4';
    load('data_cos4_eps10_VACAg.mat');   
%     load('data_cos4_eps20_VACAg_Nth128.mat');   
%     load('data_cos4_eps10_WATERAg.mat');   
%     load('data_cos4_eps20_WATERAg_Nth128.mat');   
end



fprintf('-------------\n');
fprintf('a = %g  b = %g  Eps_max = %g\n',a,b,Eps_max);
fprintf('N_theta = %d N = %d M = %d N_eps = %d\n',N_theta,N,M,N_eps);
fprintf('Material: outer = %s, inner = %s\n',OUT,IN);
fprintf('\n');


figure(1);
contourf(lambda,epsvec,U_norm,40);
colorbar;
xlabel('$\lambda$','interpreter','latex');
ylabel('$\epsilon$','interpreter','latex');
title('$|U|_2$ versus $\lambda$ and $\epsilon$','interpreter','latex');

% figure(1);
% contourf(lambda,epsvec,Qu_norm,40);
% colorbar;
% xlabel('$\lambda$','interpreter','latex');
% ylabel('$\epsilon$','interpreter','latex');
% title('$|Qu|_2$ versus $\lambda$ and $\epsilon$','interpreter','latex');
% 
% figure(2);
% contourf(lambda,epsvec,Sw_norm,40);
% colorbar;
% xlabel('$\lambda$','interpreter','latex');
% ylabel('$\epsilon$','interpreter','latex');
% title('$|Sw|_2$ versus $\lambda$ and $\epsilon$','interpreter','latex');
% 
% figure(3);
% BU_norm11 = BU_norm(1,:);
% BU_ratio = BU_norm./BU_norm11;
% contourf(lambda,epsvec,BU_ratio);
% colorbar;
% xlabel('$\lambda$','interpreter','latex');
% ylabel('$\epsilon$','interpreter','latex');
% title('$|BU_{ratio}|_2$ versus $\lambda$ and $\epsilon$','interpreter','latex');
% 
% 
% if(SavePlots==1)
%     filename = sprintf('fig_IIO_Qu_%s_eps%.0flam%s%s',name,Eps_max/a*100,...
%         IN,OUT);
%     saveas(1,filename,'epsc');
%     filename = sprintf('fig_IIO_Sw_%s_eps%.0flam%s%s',name,Eps_max/a*100,...
%         IN,OUT);
%     saveas(2,filename,'epsc');
%     filename = sprintf('fig_IIO_BU_%s_eps%.0flam%s%s',name,Eps_max/a*100,...
%         IN,OUT);
%     saveas(3,filename,'epsc');
% end

saveas(1,'LSPR_2','epsc');