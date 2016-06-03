function [relerr,nplot] = compute_errors_2d_polar(nu,Gn_fe,Eps,N,N_theta)

% relerr = zeros(N+1,6);
relerr = zeros(N+1,2);
nplot = zeros(N+1,1);

% nu_oe_taylor = Gn_oe(:,0+1);
% nu_oe_pade = Gn_oe(:,0+1);
nu_fe_taylor = Gn_fe(:,0+1);
nu_fe_pade = Gn_fe(:,0+1);
% nu_tfe_taylor = Gn_tfe(:,0+1);
% nu_tfe_pade = Gn_tfe(:,0+1);

nplot(0+1) = 0;
% relerr(0+1,1) = norm(nu-nu_oe_taylor,inf)/norm(nu,inf);
% relerr(0+1,2) = norm(nu-nu_oe_pade,inf)/norm(nu,inf);
relerr(0+1,1) = norm(nu-nu_fe_taylor,inf)/norm(nu,inf);
relerr(0+1,2) = norm(nu-nu_fe_pade,inf)/norm(nu,inf);
% relerr(0+1,5) = norm(nu-nu_tfe_taylor,inf)/norm(nu,inf);
% relerr(0+1,6) = norm(nu-nu_tfe_pade,inf)/norm(nu,inf);

for n=1:N
  M = floor(n/2);
  for j=1:N_theta
%     coeff = Gn_oe(j,:).';
%     nu_oe_taylor(j) = taylorsum(coeff,Eps,n);
%     nu_oe_pade(j) = padesum(coeff,Eps,M);
    coeff = Gn_fe(j,:).';
    nu_fe_taylor(j) = taylorsum(coeff,Eps,n);
    nu_fe_pade(j) = padesum(coeff,Eps,M);
%     coeff = Gn_tfe(j,:).';
%     nu_tfe_taylor(j) = taylorsum(coeff,Eps,n);
%     nu_tfe_pade(j) = padesum(coeff,Eps,M);
  end
  nplot(n+1) = n;
%   relerr(n+1,1) = norm(nu-nu_oe_taylor,inf)/norm(nu,inf);
%   relerr(n+1,2) = norm(nu-nu_oe_pade,inf)/norm(nu,inf);
  relerr(n+1,1) = norm(nu-nu_fe_taylor,inf)/norm(nu,inf);
  relerr(n+1,2) = norm(nu-nu_fe_pade,inf)/norm(nu,inf);
%   relerr(n+1,5) = norm(nu-nu_tfe_taylor,inf)/norm(nu,inf);
%   relerr(n+1,6) = norm(nu-nu_tfe_pade,inf)/norm(nu,inf);
end

fprintf('\n');
% fprintf('n  OE(T)  OE(P)  FE(T)  FE(P)  TFE(T)  TFE(P)\n');
% fprintf('---------------------------------------------\n');
% for n=0:N
%   fprintf('%d  %g  %g  %g  %g  %g  %g\n',nplot(n+1),...
%       relerr(n+1,1),relerr(n+1,2),relerr(n+1,3),relerr(n+1,4),...
%       relerr(n+1,5),relerr(n+1,6));
% end

fprintf('n  FE(T)  FE(P)  \n');
fprintf('-----------------\n');
for n=0:N
  fprintf('%d  %g  %g\n',nplot(n+1), relerr(n+1,1),relerr(n+1,2));
end
return;