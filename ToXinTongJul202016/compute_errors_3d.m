function [relerr,nplot] = compute_errors_3d(nu,Gn_oe,Gn_fe,Gn_tfe,Eps,N,Nx1,Nx2)

relerr = zeros(N+1,6);
nplot = zeros(N+1,1);

nu_oe_taylor = Gn_oe(:,:,0+1);
nu_oe_pade = Gn_oe(:,:,0+1);
nu_fe_taylor = Gn_fe(:,:,0+1);
nu_fe_pade = Gn_fe(:,:,0+1);
nu_tfe_taylor = Gn_tfe(:,:,0+1);
nu_tfe_pade = Gn_tfe(:,:,0+1);

nplot(0+1) = 0;
relerr(0+1,1) = norm(nu-nu_oe_taylor,inf)/norm(nu,inf);
relerr(0+1,2) = norm(nu-nu_oe_pade,inf)/norm(nu,inf);
relerr(0+1,3) = norm(nu-nu_fe_taylor,inf)/norm(nu,inf);
relerr(0+1,4) = norm(nu-nu_fe_pade,inf)/norm(nu,inf);
relerr(0+1,5) = norm(nu-nu_tfe_taylor,inf)/norm(nu,inf);
relerr(0+1,6) = norm(nu-nu_tfe_pade,inf)/norm(nu,inf);

for n=1:N
  M = floor(n/2);
  for j1=1:Nx1
    for j2=1:Nx2
      coeff = zeros(N+1,1);
      for ell=0:N
        coeff(ell+1) = Gn_oe(j1,j2,ell+1);
      end
      nu_oe_taylor(j1,j2) = taylorsum(coeff,Eps,n);
      nu_oe_pade(j1,j2) = padesum(coeff,Eps,M);
      for ell=0:N
        coeff(ell+1) = Gn_fe(j1,j2,ell+1);
      end
      nu_fe_taylor(j1,j2) = taylorsum(coeff,Eps,n);
      nu_fe_pade(j1,j2) = padesum(coeff,Eps,M);
      for ell=0:N
        coeff(ell+1) = Gn_tfe(j1,j2,ell+1);
      end
      nu_tfe_taylor(j1,j2) = taylorsum(coeff,Eps,n);
      nu_tfe_pade(j1,j2) = padesum(coeff,Eps,M);
    end
  end
  nplot(n+1) = n;
  relerr(n+1,1) = norm(nu-nu_oe_taylor,inf)/norm(nu,inf);
  relerr(n+1,2) = norm(nu-nu_oe_pade,inf)/norm(nu,inf);
  relerr(n+1,3) = norm(nu-nu_fe_taylor,inf)/norm(nu,inf);
  relerr(n+1,4) = norm(nu-nu_fe_pade,inf)/norm(nu,inf);
  relerr(n+1,5) = norm(nu-nu_tfe_taylor,inf)/norm(nu,inf);
  relerr(n+1,6) = norm(nu-nu_tfe_pade,inf)/norm(nu,inf);
end

return;