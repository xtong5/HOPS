% plot
clear all

SavePlots=0;


load('data_expcos_eps20_VACVAC_1to10.mat');
m = M;
lambda_all(1:m) = lambda;
U_norm1(1:m) = U_norm(end,:); W_norm1(1:m) = W_norm(end,:);
BU_norm1(1:m) = BU_norm(end,:); BU_norm11(1:m) = BU_norm(1,:);
Gn_U_norm1(1:m) = Gn_U_norm(end,:); Gn_W_norm1(1:m) = Gn_W_norm(end,:);

load('data_expcos_eps20_VACVAC_10to20.mat');
lambda_all(m+1:m+M) = lambda;
U_norm1(m+1:m+M) = U_norm(end,:); W_norm1(m+1:m+M) = W_norm(end,:);
BU_norm1(m+1:m+M) = BU_norm(end,:); BU_norm11(m+1:m+M) = BU_norm(1,:);
Gn_U_norm1(m+1:m+M) = Gn_U_norm(end,:); Gn_W_norm1(m+1:m+M) = Gn_W_norm(end,:);
m = m+M;

load('data_expcos_eps20_VACVAC_20to30.mat');
lambda_all(m+1:m+M) = lambda;
U_norm1(m+1:m+M) = U_norm(end,:); W_norm1(m+1:m+M) = W_norm(end,:);
BU_norm1(m+1:m+M) = BU_norm(end,:); BU_norm11(m+1:m+M) = BU_norm(1,:);
Gn_U_norm1(m+1:m+M) = Gn_U_norm(end,:); Gn_W_norm1(m+1:m+M) = Gn_W_norm(end,:);
m = m+M;

load('data_expcos_eps20_VACVAC_30to40.mat');
lambda_all(m+1:m+M) = lambda;
U_norm1(m+1:m+M) = U_norm(end,:); W_norm1(m+1:m+M) = W_norm(end,:);
BU_norm1(m+1:m+M) = BU_norm(end,:); BU_norm11(m+1:m+M) = BU_norm(1,:);
Gn_U_norm1(m+1:m+M) = Gn_U_norm(end,:); Gn_W_norm1(m+1:m+M) = Gn_W_norm(end,:);
m = m+M;

load('data_expcos_eps20_VACVAC_40to50.mat');
lambda_all(m+1:m+M) = lambda;
U_norm1(m+1:m+M) = U_norm(end,:); W_norm1(m+1:m+M) = W_norm(end,:);
BU_norm1(m+1:m+M) = BU_norm(end,:); BU_norm11(m+1:m+M) = BU_norm(1,:);
Gn_U_norm1(m+1:m+M) = Gn_U_norm(end,:); Gn_W_norm1(m+1:m+M) = Gn_W_norm(end,:);
m = m+M;

fprintf('-------------\n');
fprintf('a = %g  b = %g  Eps_max = %g\n',a,b,Eps_max);
fprintf('lambda_st = %g  lambda_end = %g\n',lambda_all(1),lambda_all(end));
fprintf('N_theta = %d N = %d M = %d N_eps = %d\n',N_theta,N,m,N_eps);
fprintf('Material: outer = %s, inner = %s\n',OUT,IN);
fprintf('\n');


BU_ratio = BU_norm1./BU_norm11;



figure(1);
norm_max = max([max(U_norm1),max(W_norm1)]);
plot(lambda_all,U_norm1,'b-o',lambda_all,W_norm1,'g-*');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|U|_2$ and $|W|_2$','interpreter','latex');
title('$|U|_2$ and $|W|_2$ versus $\lambda$','interpreter','latex');
ll = legend('$|U|_2$','$|W|_2$');
set(ll,'FontSize',16,'interpreter','latex');

figure(2);
norm_max = max([max(Gn_U_norm1),max(Gn_W_norm1)]);
plot(lambda_all,Gn_U_norm1,'b-o',lambda_all,Gn_W_norm1,'g-*');
xlabel('$\lambda$','interpreter','latex');
ylabel('$|\widetilde{U}|_2$ and $|\widetilde{W}|_2$','interpreter','latex');
title('$|\widetilde{U}|_2$ and $|\widetilde{W}|_2$ versus $\lambda$','interpreter','latex');
ll = legend('$|\widetilde{U}|_2$','$|\widetilde{W}|_2$');
set(ll,'FontSize',16,'interpreter','latex');

figure(3);
plot(lambda_all,BU_ratio,'b-o');
xlabel('$\lambda$','interpreter','latex');
ylabel('ratio');
title('$|BU|_2$ ratio versus $\lambda$','interpreter','latex');

if(SavePlots==1)
    filename = sprintf('fig_UWlam_%s_eps%.0f%s%s',name,Eps_max/a*100,...
        IN,OUT);
    saveas(1,filename,'epsc');
%     filename = sprintf('fig_index_%s_eps%.0f%s%s',name,Eps_max/a*100,...
%         IN,OUT);
%     saveas(2,filename,'epsc');
    filename = sprintf('fig_BUratiolam_%s_eps%.0f%s%s',name,Eps_max/a*100,...
        IN,OUT);
    saveas(3,filename,'epsc');
end