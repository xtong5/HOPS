clear all; close all;

PlotNum = 3;

plot_gbar = 0;
plot_ab = 0;

N = 100;
%Eps = 0.3;
Eps = 0.2;

gbar = 1.0;
b = 2.0;
c = 0.5;

dth = 2*pi/N;

th = [0:dth:2*pi];

if(PlotNum==1)
  rg = gbar + Eps*exp(cos(th));
elseif(PlotNum==2)
  rg = gbar + Eps*cos(2*th);
elseif(PlotNum==3)
  rg = gbar + Eps*cos(4*th);
else
  rg = gbar + 0*th;
end

xxgbar = gbar.*cos(th);
yygbar = gbar.*sin(th);

xxg = rg.*cos(th);
yyg = rg.*sin(th);

xxb = b.*cos(th);
yyb = b.*sin(th);

xxc = c.*cos(th);
yyc = c.*sin(th);


L = b + Eps;

plot(xxg,yyg,'k-');
hold on;
%fill(xxa,yya,'blue');
fill(xxg,yyg,[192.0/256.0 192.0/256.0 192.0/256.0]);
hold on;
if(plot_gbar==1)
  plot(xxgbar,yygbar,'b-.');
  hold on;
end


set(gca,'xtick',[]);
set(gca,'ytick',[]);

aa = 1.4;
M = 0.4;

% WaveColor = 'Green';
% 
% a = aa - 0.1; x = [-a-M;-a+M]; y = [a-M;a+M];
%   line(x,y,'LineWidth',2,'Color',WaveColor);
% a = aa; x = [-a-M;-a+M]; y = [a-M;a+M];
%   line(x,y,'LineWidth',2,'Color',WaveColor);
% a = aa + 0.1; x = [-a-M;-a+M]; y = [a-M;a+M];
%   line(x,y,'LineWidth',2,'Color',WaveColor);
% 
% hold on;
% a = aa; start = [-a-M;a+M]; stop = [-a+M;a-M];
%   arrow(start,stop,6.0,WaveColor);
% hold on;

% text(stop(1),start(2)+.2,'$(\alpha,-\gamma^u)$',...
%     'interpreter','latex','FontSize',24);
% 
axis equal;
axis([-L L -L L]);

% text(-0.2,0,'$S^w$','interpreter','latex','FontSize',40);
text(-1.2,-1.2,'$S^u$','interpreter','latex','FontSize',40);
text(-1.2,0,'$S^w$','interpreter','latex','FontSize',40);


% text(0.8,-0.6,'$\leftarrow r = \bar{g}$','Color','blue','interpreter','latex','FontSize',20);
% text(0.25,1.1,'$\leftarrow r = \bar{g} +g$','interpreter','latex','FontSize',20);

text(1.5,1.7,'$S^o$','color','r','interpreter','latex','FontSize',40);
text(-0.2,0,'$S_i$','color','m','interpreter','latex','FontSize',40);

if(plot_ab==1)
  plot(xxb,yyb,'r--');
  x1 = (b+0.2)*cos(pi/4);
  y1 = (b+0.2)*sin(pi/4);
  text(x1,y1,'$\mathcal{B}$','interpreter','latex','FontSize',40);
  x2 = (b+0.5)*cos(5*pi/4);
  y2 = (b+0.5)*sin(5*pi/4);
  text(x2,y2,'$\Omega$','interpreter','latex','FontSize',40);
end

plot(xxb,yyb,'r--');
text(-0.4,1.8,'$r=R^o$','Color','r','interpreter','latex','FontSize',24);

plot(xxc,yyc,'m--');
text(-0.4,0.6,'$r=R_i$','Color','m','interpreter','latex','FontSize',24);
saveas(gca,'boplot0_bdry2','epsc');

% saveas(gca,'boplot0','epsc');
% saveas(gca,'boplot0_1','epsc');
% saveas(gca,'boplot0_bdry1','epsc');