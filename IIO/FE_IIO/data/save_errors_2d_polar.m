function [relerr_taylor,relerr_pade,nplot] = save_errors_2d_polar(nu,Gn_fe,Eps,N,N_theta)

relerr_taylor = zeros(N+1,1);
relerr_pade = zeros(N+1,1);
nplot = zeros(N+1,1);

nu_fe_taylor = Gn_fe(:,0+1);
nu_fe_pade = Gn_fe(:,0+1);

nplot(0+1) = 0;
relerr_taylor(1) = norm(nu-nu_fe_taylor,inf)/norm(nu,inf);
relerr_pade(1) = norm(nu-nu_fe_pade,inf)/norm(nu,inf);


for n=1:N
  M = floor(n/2);
  for j=1:N_theta
    coeff = Gn_fe(j,:).';
    nu_fe_taylor(j) = taylorsum(coeff,Eps,n);
    nu_fe_pade(j) = padesum(coeff,Eps,M);
  end
  nplot(n+1) = n;
  relerr_taylor(n+1) = norm(nu-nu_fe_taylor,inf)/norm(nu,inf);
  relerr_pade(n+1)= norm(nu-nu_fe_pade,inf)/norm(nu,inf);
end


