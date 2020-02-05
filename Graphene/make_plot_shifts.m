
close all

SavePlots = 1;

% load('data_IIO_shifts_radius_VACVAC_0.025to1_ALMA32.mat');
% load('data_IIO_shifts_radius_VACVAC_0.025to1_ALMA.mat');
% load('data_IIO_shifts_Ngbars10_Nlamb51_ALMA.mat');
% load('data_IIO_shifts_Ngbars5_Nlamb101_ALMAhalf2.mat');
load('data_IIO_shifts_Nlamb101_ALMA.mat');
% load('data_IIO_shifts_Nlamb101_BFPV.mat');
% fprintf('shift of Qu = %g, shift of Sw = %g\n', d_shift_Qu(end),d_shift_Sw(end));

figure(1)
plot(g_bar,d_shift_Qu,'b-o',g_bar,d_shift_Sw,'r-*','LineWidth',2);
xlabel('$\bar{g}$','interpreter','latex','FontSize',20);
ylabel('Shift','FontSize',20);
title('Shift versus $\bar{g}$','interpreter','latex','FontSize',20);
ll = legend('Shift Qu','Shift Sw');
set(ll,'FontSize',20,'interpreter','latex','Location','best');


if(SavePlots==1)
    filename = 'Fig_shiftVSradius_VACVAC_ALMA';
%     filename = 'Fig_shiftVSradius_VACVAC_BFPV';
    saveas(1,filename,'epsc');
end


% load('data_cos4_eps20_VACVAC_34to40_gbar1_BFPV.mat');
% load('data_IIO_cos4_eps20_VACVAC_34to39_gbar1_BFPV.mat');

% load('data_IIO_cos4_eps20_VACVAC_35to40_gbar01_ALMA.mat');

% 
% [Qu_flat,Qu_index_flat] = max(Qu_norm(1,:));[Sw_flat,Sw_index_flat] = max(Sw_norm(1,:));
% [Qu_shift,Qu_index_shift] = max(Qu_norm(end,:));[Sw_shift,Sw_index_shift] = max(Sw_norm(end,:));
% d_shift_Qu = abs(lambda(Qu_index_flat) - lambda(Qu_index_shift));
% d_shift_Sw = abs(lambda(Sw_index_flat) - lambda(Sw_index_shift));
% fprintf('shift of Qu = %g, shift of Sw = %g\n', d_shift_Qu,d_shift_Sw);
