% Nplot = [0 4 8 12 16];
SavePlots=0;
Nplot = [0 2 4 6 8 10 12 14 16];
N = size(Nplot,2);
mode = 2; %1=U, 2=IIOU/DNOU, 3=W, 4=IIOW/DNOW
% mode_op = 1; Operater = 'IIO_fe'; %1=IIO, 2=DNO
mode_op = 3; Operater = 'IIO_new_fe';
% mode_op = 2; Operater = 'DNO_fe'; 
% mode_sing = ''; %empty, _sing12, _sing16
% mode_sing = '_sing12'; 
mode_sing = '_sing16'; 

name1 = sprintf('errors_%s_Eps_0.005%s.mat',Operater,mode_sing);
name2 = sprintf('errors_%s_Eps_0.01%s.mat',Operater,mode_sing);
name3 = sprintf('errors_%s_Eps_0.05%s.mat',Operater,mode_sing);
name4 = sprintf('errors_%s_Eps_0.1%s.mat',Operater,mode_sing);
    
if mode == 1
    load(name1);relerr1 = relerr_padeU;
    load(name2);relerr2 = relerr_padeU;
    load(name3);relerr3 = relerr_padeU;
    load(name4);relerr4 = relerr_padeU;
end

if mode == 2
    if mode_op == 2
        load(name1);relerr1 = relerr_padeDNOU;
        load(name2);relerr2 = relerr_padeDNOU;
        load(name3);relerr3 = relerr_padeDNOU;
        load(name4);relerr4 = relerr_padeDNOU;
    else
    load(name1);relerr1 = relerr_padeIIOU;
    load(name2);relerr2 = relerr_padeIIOU;
    load(name3);relerr3 = relerr_padeIIOU;
    load(name4);relerr4 = relerr_padeIIOU;
    end
end


if mode == 3
    load(name1);relerr1 = relerr_padeW;
    load(name2);relerr2 = relerr_padeW;
    load(name3);relerr3 = relerr_padeW;
    load(name4);relerr4 = relerr_padeW;
end

if mode == 4
    if mode_op == 2
        load(name1);relerr1 = relerr_padeDNOW;
        load(name2);relerr2 = relerr_padeDNOW;
        load(name3);relerr3 = relerr_padeDNOW;
        load(name4);relerr4 = relerr_padeDNOW;
    else
    load(name1);relerr1 = relerr_padeIIOW;
    load(name2);relerr2 = relerr_padeIIOW;
    load(name3);relerr3 = relerr_padeIIOW;
    load(name4);relerr4 = relerr_padeIIOW;
    end
end

for i=1:N
    error1(i)=relerr1(Nplot(i)+1);
    error2(i)=relerr2(Nplot(i)+1);
    error3(i)=relerr3(Nplot(i)+1);
    error4(i)=relerr4(Nplot(i)+1);
end

figure(1);
semilogy(Nplot,error1,'c-d',Nplot,error2,'y-+',Nplot,error3,'r-o',Nplot,error4,'b-*');
xlabel('$N$','interpreter','latex');
ylabel('Relative Error','interpreter','latex');
title('Relative Error versus $N$','interpreter','latex');
ll = legend('$\varepsilon=0.005$','$\varepsilon=0.01$','$\varepsilon=0.05$','$\varepsilon=0.1$');
set(ll,'FontSize',10,'interpreter','latex');

eps = [0.005 0.01 0.05 0.1];
Error0(1) = relerr1(1);Error0(2) = relerr2(1);Error0(3) = relerr3(1);Error0(4) = relerr4(1);
Error4(1) = relerr1(1+4);Error4(2) = relerr2(1+4);Error4(3) = relerr3(1+4);Error4(4) = relerr4(1+4);
Error8(1) = relerr1(1+8);Error8(2) = relerr2(1+8);Error8(3) = relerr3(1+8);Error8(4) = relerr4(1+8);
Error12(1) = relerr1(1+12);Error12(2) = relerr2(1+12);Error12(3) = relerr3(1+12);Error12(4) = relerr4(1+12);
Error16(1) = relerr1(1+16);Error16(2) = relerr2(1+16);Error16(3) = relerr3(1+16);Error16(4) = relerr4(1+16);

figure(2);
% semilogy(eps,Error0,'c-d',eps,Error4,'y-+',eps,Error8,'r-o',eps,Error12,'b-*',eps,Error16,'k-^');
loglog(eps,Error0,'c-d',eps,Error4,'y-+',eps,Error8,'r-o',eps,Error12,'b-*',eps,Error16,'k-^');
xlabel('$\varepsilon$','interpreter','latex');
ylabel('Relative Error','interpreter','latex');
title('Relative Error versus $\varepsilon$','interpreter','latex');
ll = legend('$N=0$','$N=4$','$N=8$','$N=12$','$N=16$');
set(ll,'FontSize',10,'interpreter','latex');

if(SavePlots==1)
    filename = sprintf('fig_ErrrorVsN_%s%s', Operater,mode_sing);
    saveas(1,filename,'epsc');
    filename = sprintf('fig_ErrrorVsEps_%s%s', Operater,mode_sing);
    saveas(2,filename,'epsc');
end