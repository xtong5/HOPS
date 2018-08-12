function [err_taylor,err_pade,nplot] = save_abs_errors_2d_polar(nu,Gn_tfe,Eps,N,N_theta)

err_taylor = zeros(N+1,1);
err_pade = zeros(N+1,1);
nplot = zeros(N+1,1);

nu_tfe_taylor = Gn_tfe(:,0+1);
nu_tfe_pade = Gn_tfe(:,0+1);

nplot(0+1) = 0;
err_taylor(1) = norm(nu-nu_tfe_taylor,inf);
err_pade(1) = norm(nu-nu_tfe_pade,inf);


for n=1:N
  M = floor(n/2);
  for j=1:N_theta
    coeff = Gn_tfe(j,:).';
    nu_tfe_taylor(j) = taylorsum(coeff,Eps,n);
    nu_tfe_pade(j) = padesum(coeff,Eps,M);
  end
  nplot(n+1) = n;
  err_taylor(n+1) = norm(nu-nu_tfe_taylor,inf);
  err_pade(n+1)= norm(nu-nu_tfe_pade,inf);
end


