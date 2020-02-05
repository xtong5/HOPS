function [] = save_plots(SavePlots,nplot,relerr)

%
% Taylor Summation
%

% FE

figure(11);
clf;
hh11 = gca;
semilogy(nplot,relerr(:,1),'k-^',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('FE(Taylor)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh11,'fontsize',16);
% if(SavePlots==1)
%   saveas(hh11,'err_fe_taylor','eps');
% end

% TFE

figure(13);
clf;
hh13 = gca;
semilogy(nplot,relerr(:,3),'k-o',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('TFE(Taylor)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh13,'fontsize',16);
% if(SavePlots==1)
%   saveas(hh13,'err_tfe_taylor','eps');
% end

%
% Pade Summation
%

% FE

figure(21);
clf;
hh21 = gca;
semilogy(nplot,relerr(:,2),'k-^',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('FE(Pade)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh21,'fontsize',16);
% if(SavePlots==1)
%   saveas(hh21,'err_fe_pade','eps');
% end

% TFE

figure(23);
clf;
hh23 = gca;
semilogy(nplot,relerr(:,4),'k-o',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('TFE(Pade)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh23,'fontsize',16);
% if(SavePlots==1)
%   saveas(hh23,'err_tfe_pade','eps');
% end

%
% Taylor and Pade
%

% FE

figure(101);
clf;
hh101 = gca;
semilogy(nplot,relerr(:,1),'k-^',nplot,relerr(:,2),'k-^',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('FE(Taylor)','FE(Pade)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh101,'fontsize',16);
% if(SavePlots==1)
%   saveas(hh101,'err_fe_taylorpade','eps');
% end

% TFE

figure(103);
clf;
hh103 = gca;
semilogy(nplot,relerr(:,3),'k-o',nplot,relerr(:,4),'k-o',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('TFE(Taylor)','TFE(Pade)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh103,'fontsize',16);
% if(SavePlots==1)
%   saveas(hh103,'err_tfe_taylorpade','eps');
% end


%% All
figure(1234);
clf;
hh1234 = gca;
semilogy(nplot,relerr(:,1),'k-v',nplot,relerr(:,2),'k-^',...
    nplot,relerr(:,3),'k-s',nplot,relerr(:,4),'k-d',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('FE(Taylor)','FE(Pade)','TFE(Taylor)','TFE(Pade)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh1234,'fontsize',16);





