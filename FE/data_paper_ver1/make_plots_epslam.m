% plot
clear all
close all
SavePlots = 1;

mode = 4; % choose functions

if mode == 1
    name = 'expcos';
%     load('data_expcos_eps10_VACAg.mat');
%     load('data_expcos_eps10_WATERAg.mat');
%     load('data_expcos_eps20_VACAg.mat');
    load('data_expcos_eps20_WATERAg.mat');
end
if mode == 2
    name = 'cos2';
%     load('data_cos2_eps20_VACAg.mat');
%     load('data_cos2_eps20_WATERAg.mat');
%     load('data_cos2_eps10_VACAg.mat');
    load('data_cos2_eps10_WATERAg.mat');
end
if mode == 4
    name = 'cos4';
%     load('data_cos4_eps10_VACAg.mat');
%     load('data_cos4_eps10_WATERAg.mat');
%     load('data_cos4_eps20_VACAgNt128N16.mat');
    load('data_cos4_eps20_WATERAgNt128N16.mat');
end

fprintf('-------------\n');
fprintf('a = %g  Eps_max = %g\n',a,Eps_max);
fprintf('N_theta = %d N = %d N_lamb = %d N_eps = %d\n',N_theta,N,M,N_eps);
fprintf('Material: outer = %s, inner = %s\n',OUT,IN);
fprintf('\n');

figure(1);
hh = gca;
contourf(lambda,epsvec,Gn_U_norm,40);
colorbar;
xlabel('$\lambda$','interpreter','latex');
ylabel('$\varepsilon$','interpreter','latex');
title('$|\widetilde{U}|_2$ versus $\lambda$ and $\varepsilon$','interpreter','latex');
set(hh,'FontSize',16);

figure(2);
hh = gca;
contourf(lambda,epsvec,Gn_W_norm,40);
colorbar;
xlabel('$\lambda$','interpreter','latex');
ylabel('$\varepsilon$','interpreter','latex');
title('$|\widetilde{W}|_2$ versus $\lambda$ and $\varepsilon$','interpreter','latex');
set(hh,'FontSize',16);

if(SavePlots==1)
    filename = sprintf('fig_U_%s_eps%.0flam%s%s',name,Eps_max/a*100,...
        IN,OUT);
    saveas(1,filename,'epsc');
    filename = sprintf('fig_W_%s_eps%.0flam%s%s',name,Eps_max/a*100,...
        IN,OUT);
    saveas(2,filename,'epsc');
end