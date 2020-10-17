% plot
clear all

load('data_IIO_cos4_eps20_VACVAC_35to40_gbar01_ALMA.mat');
Iu_norm1 = Iu_norm(end,:); Iw_norm1 = Iw_norm(end,:); % last slice or epsilon=max
Qu_norm1 = Qu_norm(end,:); Sw_norm1 = Sw_norm(end,:);
[maxQ1,indQ1]=max(Qu_norm1);[maxS1,indS1]=max(Sw_norm1); % find lspr 
lamb_Q1=lambda(indQ1);lamb_S1=lambda(indS1); % find index of lambda

load('data_nIIO_cos4_eps20_VACVAC_35to40_gbar01_ALMA.mat');
Iu_norm2 = Iu_norm(end,:); Iw_norm2 = Iw_norm(end,:); 
Qu_norm2 = Qu_norm(end,:); Sw_norm2 = Sw_norm(end,:);
[maxQ2,indQ2]=max(Qu_norm2);[maxS2,indS2]=max(Sw_norm2);
lamb_Q2=lambda(indQ2);lamb_S2=lambda(indS2);


% Qu_norm0 = Qu_norm(1,:); Sw_norm0 = Sw_norm(1,:);

% Iu_norm1-Iu_norm2
% Iw_norm1-Iw_norm2
% Qu_norm1-Qu_norm2

fprintf('maxQ1 = %g  indQ1 = %g  lambdaQ1 = %g\n',maxQ1,indQ1,lamb_Q1);
fprintf('maxS1 = %g  indS1 = %g  lambdaS1 = %g\n',maxS1,indS1,lamb_S1);
fprintf('maxQ2 = %g  indQ2 = %g  lambdaQ2 = %g\n',maxQ2,indQ2,lamb_Q2);
fprintf('maxS2 = %g  indS2 = %g  lambdaS2 = %g\n',maxS2,indS2,lamb_S2);

