function [] = make_plots(SavePlots,nplot,relerr)

%
% Taylor Summation
%

% OE

figure(11);
clf;
hh11 = gca;
semilogy(nplot,relerr(:,1),'k-^',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('OE(Taylor)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh11,'fontsize',16);
if(SavePlots==1)
  saveas(hh11,'err_oe_taylor','eps');
end

% FE

figure(12);
clf;
hh12 = gca;
semilogy(nplot,relerr(:,3),'k-*',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('FE(Taylor)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh12,'fontsize',16);
if(SavePlots==1)
  saveas(hh12,'err_fe_taylor','eps');
end

% TFE

figure(13);
clf;
hh13 = gca;
semilogy(nplot,relerr(:,5),'k-o',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('TFE(Taylor)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh13,'fontsize',16);
if(SavePlots==1)
  saveas(hh13,'err_tfe_taylor','eps');
end

% OE and FE

figure(14);
clf;
hh14 = gca;
semilogy(nplot,relerr(:,1),'k-^',...
    nplot,relerr(:,3),'k-*',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('OE(Taylor)','FE(Taylor)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh14,'fontsize',16);
if(SavePlots==1)
  saveas(hh14,'err_oefe_taylor','eps');
end

% FE, OE, and TFE

figure(15);
clf;
hh15 = gca;
semilogy(nplot,relerr(:,1),'k-^',...
    nplot,relerr(:,3),'k-*',...
    nplot,relerr(:,5),'k-o',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('OE(Taylor)','FE(Taylor)','TFE(Taylor)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh15,'fontsize',16);
if(SavePlots==1)
  saveas(hh15,'err_oefetfe_taylor','eps');
end

%
% Pade Summation
%

% OE

figure(21);
clf;
hh21 = gca;
semilogy(nplot,relerr(:,2),'k-^',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('OE(Pade)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh21,'fontsize',16);
if(SavePlots==1)
  saveas(hh21,'err_oe_pade','eps');
end

% FE

figure(22);
clf;
hh22 = gca;
semilogy(nplot,relerr(:,4),'k-*',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('FE(Pade)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh22,'fontsize',16);
if(SavePlots==1)
  saveas(hh22,'err_fe_pade','eps');
end

% TFE

figure(23);
clf;
hh23 = gca;
semilogy(nplot,relerr(:,6),'k-o',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('TFE(Pade)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh23,'fontsize',16);
if(SavePlots==1)
  saveas(hh23,'err_tfe_pade','eps');
end

% OE and FE

figure(24);
clf;
hh24 = gca;
semilogy(nplot,relerr(:,2),'k-^',...
    nplot,relerr(:,4),'k-*',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('OE(Pade)','FE(Pade)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh24,'fontsize',16);
if(SavePlots==1)
  saveas(hh23,'err_oefe_pade','eps');
end

% OE, FE, and TFE

figure(25);
clf;
hh25 = gca;
semilogy(nplot,relerr(:,2),'k-^',...
    nplot,relerr(:,4),'k-*',...
    nplot,relerr(:,6),'k-o',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('OE(Pade)','FE(Pade)','TFE(Pade)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh25,'fontsize',16);
if(SavePlots==1)
  saveas(hh25,'err_oefetfe_pade','eps');
end

%
% Taylor and Pade
%

% OE

figure(101);
clf;
hh101 = gca;
semilogy(nplot,relerr(:,1),'k-^',nplot,relerr(:,2),'k-^',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('OE(Taylor)','OE(Pade)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh101,'fontsize',16);
if(SavePlots==1)
  saveas(hh101,'err_oe_taylorpade','eps');
end

% FE

figure(102);
clf;
hh102 = gca;
semilogy(nplot,relerr(:,2),'k-*',nplot,relerr(:,4),'k-*',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('FE(Taylor)','FE(Pade)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh102,'fontsize',16);
if(SavePlots==1)
  saveas(hh102,'err_fe_taylorpade','eps');
end

% TFE

figure(103);
clf;
hh103 = gca;
semilogy(nplot,relerr(:,3),'k-o',nplot,relerr(:,6),'k-o',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('TFE(Taylor)','TFE(Pade)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh103,'fontsize',16);
if(SavePlots==1)
  saveas(hh103,'err_tfe_taylorpade','eps');
end

% OE and OE

figure(104);
clf;
hh104 = gca;
semilogy(nplot,relerr(:,1),'k-^',nplot,relerr(:,2),'k-^',...
    nplot,relerr(:,3),'k-*',nplot,relerr(:,4),'k-*',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('OE(Taylor)','OE(Pade)','FE(Taylor)','FE(Pade)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh104,'fontsize',16);
if(SavePlots==1)
  saveas(hh104,'err_oefe_taylorpade','eps');
end

% OE, FE, and TFE

figure(105);
clf;
hh105 = gca;
semilogy(nplot,relerr(:,1),'k-^',nplot,relerr(:,2),'k-^',...
    nplot,relerr(:,3),'k-*',nplot,relerr(:,4),'k-*',...
    nplot,relerr(:,5),'k-o',nplot,relerr(:,6),'k-o',...
    'LineWidth',2,'MarkerSize',8,'MarkerFaceColor','k');
ll = legend('OE(Taylor)','OE(Pade)','FE(Taylor)','FE(Pade)','TFE(Taylor)','TFE(Pade)');
set(ll,'interpreter','latex');
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(hh105,'fontsize',16);
if(SavePlots==1)
  saveas(hh105,'err_oefetfe_taylorpade','eps');
end