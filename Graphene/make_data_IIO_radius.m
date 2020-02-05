% MAKE_DATA_IIO.m 
%
% compute and save the date to file 
%
% XT 1/20

clear all
warning('off')
save_data = 1;

L = 2*pi;
N = 16;
N_theta = 64; N_eps = 2;
theta = (L/N_theta)*[0:N_theta-1]';
Y_p = 1i*3.4.*ones(N_theta,1); Z_p = -1i*3.4.*ones(N_theta,1);
N_lambda = 101; lambda_low = 32; lambda_high = 45;

OUT = 'VACUUM'; IN = 'VACUUM'; mu = 0.4;

% f0 = exp(cos(theta)); name = 'expcos';
% f2 = cos(2*theta); name = 'cos2';    
f4 = cos(4*theta); name = 'cos4';
 
g_bar = [0.025, 0.05, 0.1, 0.2, 0.5, 1];
N_gbar = size(g_bar,2);

d_shift_Qu = zeros(N_gbar,1);
d_shift_Sw = zeros(N_gbar,1);
Qu_norm_flat = zeros(N_lambda,N_gbar); Sw_norm_flat = zeros(N_lambda,N_gbar);
Qu_norm_shift = zeros(N_lambda,N_gbar); Sw_norm_shift = zeros(N_lambda,N_gbar);
lambda = linspace(lambda_low,lambda_high,N_lambda);

tic
for i = 1:N_gbar
    gbar = g_bar(i); b = 10.*gbar; Eps_max = 0.2.*gbar;
%     [Qu_norm,Sw_norm,BU_norm] = FE_app_IIO_ALMA_nosave(f4,N_theta,theta,...
%     gbar,b,Y_p,Z_p, N,Eps_max,N_eps,OUT,IN,mu,lambda_low,lambda_high,N_lambda);
    [Qu_norm,Sw_norm,BU_norm] = FE_app_IIO_BFPV_nosave(f4,N_theta,theta,...
    gbar,b,Y_p,Z_p, N,Eps_max,N_eps,OUT,IN,mu,lambda_low,lambda_high,N_lambda);
    Qu_norm_flat(:,i) = Qu_norm(1,:).'; Sw_norm_flat(:,i) = Sw_norm(1,:).';
    Qu_norm_shift(:,i) = Qu_norm(end,:).'; Sw_norm_shift(:,i) = Sw_norm(end,:).';
    [Qu_flat,Qu_index_flat] = max(Qu_norm(1,:));[Sw_flat,Sw_index_flat] = max(Sw_norm(1,:));
    [Qu_shift,Qu_index_shift] = max(Qu_norm(end,:));[Sw_shift,Sw_index_shift] = max(Sw_norm(end,:));
    d_shift_Qu(i) = abs(lambda(Qu_index_flat) - lambda(Qu_index_shift));
    d_shift_Sw(i) = abs(lambda(Sw_index_flat) - lambda(Sw_index_shift));
end
toc

filename = sprintf('data_IIO_shifts_Nlamb%g_BFPV',N_lambda);
save(filename,'Qu_norm_flat','Sw_norm_flat','Qu_norm_shift','Sw_norm_shift','d_shift_Qu','d_shift_Sw','g_bar',...
    'N_eps','N_lambda');


