function [un,wn] = field_tfe_helmholtz_twolayer(zeta_n,psi_n,f,tau2,...
    p,alphap,beta_up,beta_wp,eep,eem,Dy,a,b,Nx,Ny,N)

un = zeros(Nx,Ny+1,N+1);
wn = zeros(Nx,Ny+1,N+1);
zetahat_n = zeros(Nx,N+1);
psihat_n = zeros(Nx,N+1);

ku2 = alphap(0+1)^2 + beta_up(0+1)^2;
kw2 = alphap(0+1)^2 + beta_wp(0+1)^2;

ell_u_top = 0 + 1;
ell_u_bottom = Ny + 1;
ell_w_top = 0 + 1;
ell_w_bottom = Ny + 1;
for n=0:N
  zetahat_n(:,n+1) = fft(eem.*zeta_n(:,n+1));
  psihat_n(:,n+1) = fft(eem.*psi_n(:,n+1));
end
f_x = real(ifft( (1i*p).*fft(f) ));
ff = f.*f;
ff_x = f.*f_x;
f_xf_x = f_x.*f_x;
fff_x = f.*(f.*f_x);
ff_xf_x = f.*(f_x.*f_x);

ll = [0:Ny]';
y_u_min = 0.0; y_u_max = a;
y_u = ((y_u_max-y_u_min)/2.0)*(cos(pi*ll/Ny) - 1.0) + y_u_max;
y_w_min = -b; y_w_max = 0.0;
y_w = ((y_w_max-y_w_min)/2.0)*(cos(pi*ll/Ny) - 1.0) + y_w_max;

f_full = zeros(Nx,Ny+1);
f_x_full = zeros(Nx,Ny+1);
for ell=0:Ny
  f_full(:,ell+1) = f;
  f_x_full(:,ell+1) = f_x;
end
Uhat = zeros(Nx,Ny+1);
What = zeros(Nx,Ny+1);

a_minus_y_full = zeros(Nx,Ny+1);
y_plus_b_full = zeros(Nx,Ny+1);
for j=1:Nx
  a_minus_y_full(j,:) = a - y_u.';
  y_plus_b_full(j,:) = y_w.' + b;
end

% Order zero

Delta = tau2*1i*beta_wp + 1i*beta_up;
M0Inv_11 = tau2*(1i*beta_wp)./Delta;
M0Inv_12 = 1.0./Delta;
M0Inv_21 = -1i*beta_up./Delta;
M0Inv_22 = 1.0./Delta;
Qhat = fft( eem.*zeta_n(:,0+1) );
Rhat = fft( eem.*psi_n(:,0+1) );
uhat = M0Inv_11.*Qhat + M0Inv_12.*Rhat;
what = M0Inv_21.*Qhat + M0Inv_22.*Rhat;

for ell=0:Ny
  un(:,ell+1,0+1) = eep.*ifft( exp(1i*beta_up*y_u(ell+1)).*uhat );
  wn(:,ell+1,0+1) = eep.*ifft( exp(-1i*beta_wp*y_w(ell+1)).*what );
end

% Order n>0

for n=1:N
    
  % Form Qn
  
  Qn = psi_n(:,n+1);
  if(n>=1)
    Qn = Qn + ((a-b)/(a*b))*f.*psi_n(:,n-1+1);
  end
  if(n>=2)
    Qn = Qn - (1.0/(a*b))*ff.*psi_n(:,n-2+1);
  end
  
  % Form Fun, Jun
  
  Fun = zeros(Nx,Ny+1);
  Jun = zeros(Nx,1);
  
  A1_xx = -(2.0/a)*f_full;
  A1_xy = -(1.0/a)*(a_minus_y_full).*f_x_full;
  A1_yx = A1_xy;
  %A1_yy = 0;
  
  A2_xx = (1.0/a^2)*f_full.^2;
  A2_xy = (1.0/a^2)*(a_minus_y_full).*(f_full.*f_x_full);
  A2_yx = A2_xy;
  A2_yy = (1.0/a^2)*((a_minus_y_full).^2).*(f_x_full.^2);
  
  B1_x = -(1.0/a)*f_x_full;
  %B1_y = 0;
  
  B2_x = (1.0/a^2)*f_full.*f_x_full;
  B2_y = (1.0/a^2).*(a_minus_y_full).*(f_x_full.^2);
  
  C1 = -ku2*(2.0/a)*f_full;
  C2 = ku2*(1.0/a^2)*f_full.^2;
  
  if(n>=1)
    u_x = dxp(un(:,:,n-1+1),alphap,eep,eem,Nx,Ny);
    temp = A1_xx.*u_x;
    Fun = Fun - dxp(temp,alphap,eep,eem,Nx,Ny);
    temp = A1_yx.*u_x;
    Fun = Fun - dy(temp,Dy,a,Nx,Ny);
    temp = B1_x.*u_x;
    Fun = Fun + temp;
    
    u_y = dy(un(:,:,n-1+1),Dy,a,Nx,Ny);
    temp = A1_xy.*u_y;
    Fun = Fun - dxp(temp,alphap,eep,eem,Nx,Ny);
    %A1_yy = 0
    %B1_y = 0
    
    temp = C1.*un(:,:,n-1+1);
    Fun = Fun - temp;
    
    Su = eep.*ifft( (1i*beta_up).*fft(eem.*un(:,ell_u_top,n-1+1)) );
    Jun = Jun - (1.0/a)*f.*Su;
    
    Qn = Qn - (a/(a*b))*f.*u_y(:,ell_u_bottom);
    Qn = Qn + ((a*b)/(a*b))*f_x.*u_x(:,ell_u_bottom);
  end
  
  if(n>=2)
    u_x = dxp(un(:,:,n-2+1),alphap,eep,eem,Nx,Ny);
    temp = A2_xx.*u_x;
    Fun = Fun - dxp(temp,alphap,eep,eem,Nx,Ny);
    temp = A2_yx.*u_x;
    Fun = Fun - dy(temp,Dy,a,Nx,Ny);
    temp = B2_x.*u_x;
    Fun = Fun + temp;
    
    u_y = dy(un(:,:,n-2+1),Dy,a,Nx,Ny);
    temp = A2_xy.*u_y;
    Fun = Fun - dxp(temp,alphap,eep,eem,Nx,Ny);
    temp = A2_yy.*u_y;
    Fun = Fun - dy(temp,Dy,a,Nx,Ny);
    temp = B2_y.*u_y;
    Fun = Fun + temp;
    
    temp = C2.*un(:,:,n-2+1);
    Fun = Fun - temp;
    
    Qn = Qn + ((a-b)/(a*b))*ff_x.*u_x(:,ell_u_bottom);
    Qn = Qn - ((a*b)/(a*b))*f_xf_x.*u_y(:,ell_u_bottom);
  end
  
  if(n>=3)
    u_x = dxp(un(:,:,n-3+1),alphap,eep,eem,Nx,Ny);
    Qn = Qn - (1.0/(a*b))*fff_x.*u_x(:,ell_u_bottom);
    u_y = dy(un(:,:,n-3+1),Dy,a,Nx,Ny);
    Qn = Qn - (a/(a*b))*ff_xf_x.*u_y(:,ell_u_bottom);
  end
  
  % Form Fwn, Jwn
  
  Fwn = zeros(Nx,Ny+1);
  Jwn = zeros(Nx,1);
  
  A1_xx = (2.0/b)*f_full;
  A1_xy = -(1.0/b)*(y_plus_b_full).*f_x_full;
  A1_yx = A1_xy;
  %A1_yy = 0;
  
  A2_xx = (1.0/b^2)*f_full.^2;
  A2_xy = -(1.0/b^2)*(y_plus_b_full).*(f_full.*f_x_full);
  A2_yx = A2_xy;
  A2_yy = (1.0/b^2)*((y_plus_b_full).^2).*(f_x_full.^2);
  
  B1_x = (1.0/b)*f_x_full;
  %B1_y = 0;
  
  B2_x = (1.0/b^2)*f_full.*f_x_full;
  B2_y = -(1.0/b^2).*(y_plus_b_full).*(f_x_full.^2);
  
  C1 = kw2*(2.0/b)*f_full;
  C2 = kw2*(1.0/b^2)*f_full.^2;
  
  if(n>=1)
    w_x = dxp(wn(:,:,n-1+1),alphap,eep,eem,Nx,Ny);
    temp = A1_xx.*w_x;
    Fwn = Fwn - dxp(temp,alphap,eep,eem,Nx,Ny);
    temp = A1_yx.*w_x;
    Fwn = Fwn - dy(temp,Dy,b,Nx,Ny);
    temp = B1_x.*w_x;
    Fwn = Fwn + temp;
    
    w_y = dy(wn(:,:,n-1+1),Dy,b,Nx,Ny);
    temp = A1_xy.*w_y;
    Fwn = Fwn - dxp(temp,alphap,eep,eem,Nx,Ny);
    %A1_yy = 0
    %B1_y = 0
    
    temp = C1.*wn(:,:,n-1+1);
    Fwn = Fwn - temp;
    
    Sw = eep.*ifft( (-1i*beta_wp).*fft(eem.*wn(:,ell_w_bottom,n-1+1)) );
    Jwn = Jwn + (1.0/b)*f.*Sw;
    
    Qn = Qn - tau2*(b/(a*b))*f.*w_y(:,ell_w_top);
    Qn = Qn - tau2*((a*b)/(a*b))*f_x.*w_x(:,ell_w_top);
  end
  
  if(n>=2)
    w_x = dxp(wn(:,:,n-2+1),alphap,eep,eem,Nx,Ny);
    temp = A2_xx.*w_x;
    Fwn = Fwn - dxp(temp,alphap,eep,eem,Nx,Ny);
    temp = A2_yx.*w_x;
    Fwn = Fwn - dy(temp,Dy,b,Nx,Ny);
    temp = B2_x.*w_x;
    Fwn = Fwn + temp;
    
    w_y = dy(wn(:,:,n-2+1),Dy,b,Nx,Ny);
    temp = A2_xy.*w_y;
    Fwn = Fwn - dxp(temp,alphap,eep,eem,Nx,Ny);
    temp = A2_yy.*w_y;
    Fwn = Fwn - dy(temp,Dy,b,Nx,Ny);
    temp = B2_y.*w_y;
    Fwn = Fwn + temp;
    
    temp = C2.*wn(:,:,n-2+1);
    Fwn = Fwn - temp;
    
    Qn = Qn - tau2*((a-b)/(a*b))*ff_x.*w_x(:,ell_w_top);
    Qn = Qn + tau2*((a*b)/(a*b))*f_xf_x.*w_y(:,ell_w_top);
  end
  
  if(n>=3)
    w_x = dxp(wn(:,:,n-3+1),alphap,eep,eem,Nx,Ny);
    Qn = Qn + tau2*(1.0/(a*b))*fff_x.*w_x(:,ell_w_top);
    w_y = dy(wn(:,:,n-3+1),Dy,b,Nx,Ny);
    Qn = Qn - tau2*(b/(a*b))*ff_xf_x.*w_y(:,ell_w_top);
  end
  
  % Solve elliptic equation
  
  Funhat = zeros(Nx,Ny+1);
  Fwnhat = zeros(Nx,Ny+1);
  for ell=0:Ny
    Funhat(:,ell+1) = fft(eem.*Fun(:,ell+1));
    Fwnhat(:,ell+1) = fft(eem.*Fwn(:,ell+1));
  end
  Junhat = fft(eem.*Jun);
  Jwnhat = fft(eem.*Jwn);
  Qnhat = fft(eem.*Qn);
  
  for j=1:Nx
    Funhat_p = Funhat(j,:).';
    Fwnhat_p = Fwnhat(j,:).';
    alpha_u = 1.0;
    beta_u = 0.0;
    gamma_u = ku2 - (alphap(j))^2;
    alpha_w = 1.0;
    beta_w = 0.0;
    gamma_w = kw2 - (alphap(j))^2;
    d_b = -(-1i*beta_wp(j));
    n_b = 1.0;
    r_b = Jwnhat(j);
    dw_ab = -1.0;
    du_ab = 1.0;
    rd_ab = zetahat_n(j,n+1);
    nw_ab = -tau2;
    nu_ab = 1.0;
    rn_ab = Qnhat(j);
    d_a = -1i*beta_up(j);
    n_a = 1.0;
    r_a = Junhat(j);

    [what_p,uhat_p] = solvetwobvp_colloc(Fwnhat_p,alpha_w,beta_w,gamma_w,...
        Funhat_p,alpha_u,beta_u,gamma_u,(2.0/b)*Dy,(2.0/a)*Dy,...
        d_b,n_b,r_b,dw_ab,du_ab,nw_ab,nu_ab,rd_ab,rn_ab,d_a,n_a,r_a);
    
    Uhat(j,:) = uhat_p.';
    What(j,:) = what_p.';
  end
  
  for ell=0:Ny
    un(:,ell+1,n+1) = eep.*ifft(Uhat(:,ell+1));
    wn(:,ell+1,n+1) = eep.*ifft(What(:,ell+1));
  end

end