function [un,wn] = field_tfe_helmholtz_twolayer_3d(zeta_n,psi_n,f,tau2,...
    p1,p2,alpha1p,alpha2p,beta_up,beta_wp,eep,eem,Dy,a,b,Nx1,Nx2,Ny,N)

un = zeros(Nx1,Nx2,Ny+1,N+1);
wn = zeros(Nx1,Nx2,Ny+1,N+1);
zetahat_n = zeros(Nx1,Nx2,N+1);
psihat_n = zeros(Nx1,Nx2,N+1);

ku2 = alpha1p(0+1)^2 + alpha2p(0+1)^2 + beta_up(0+1)^2;
kw2 = alpha1p(0+1)^2 + alpha2p(0+1)^2 + beta_wp(0+1)^2;

ell_u_top = 0 + 1;
ell_u_bottom = Ny + 1;
ell_w_top = 0 + 1;
ell_w_bottom = Ny + 1;
for n=0:N
  zetahat_n(:,:,n+1) = fft2(eem.*zeta_n(:,:,n+1));
  psihat_n(:,:,n+1) = fft2(eem.*psi_n(:,:,n+1));
end
f_x1 = real(ifft2( (1i*p1).*fft2(f) ));
f_x2 = real(ifft2( (1i*p2).*fft2(f) ));
ff = f.*f;
ff_x1 = f.*f_x1;
ff_x2 = f.*f_x2;
f_x1f_x1 = f_x1.*f_x1;
f_x2f_x2 = f_x2.*f_x2;
fff_x1 = f.*(f.*f_x1);
fff_x2 = f.*(f.*f_x2);
ff_x1f_x1 = f.*(f_x1.*f_x1);
ff_x2f_x2 = f.*(f_x2.*f_x2);

ll = [0:Ny]';
y_u_min = 0.0; y_u_max = a;
y_u = ((y_u_max-y_u_min)/2.0)*(cos(pi*ll/Ny) - 1.0) + y_u_max;
y_w_min = -b; y_w_max = 0.0;
y_w = ((y_w_max-y_w_min)/2.0)*(cos(pi*ll/Ny) - 1.0) + y_w_max;

f_full = zeros(Nx1,Nx2,Ny+1);
f_x1_full = zeros(Nx1,Nx2,Ny+1);
f_x2_full = zeros(Nx1,Nx2,Ny+1);
for ell=0:Ny
  f_full(:,:,ell+1) = f;
  f_x1_full(:,:,ell+1) = f_x1;
  f_x2_full(:,:,ell+1) = f_x2;
end
Uhat = zeros(Nx1,Nx2,Ny+1);
What = zeros(Nx1,Nx2,Ny+1);

a_minus_y_full = zeros(Nx1,Nx2,Ny+1);
y_plus_b_full = zeros(Nx1,Nx2,Ny+1);
for j1=1:Nx1
  for j2=1:Nx2
    for ell=0:Ny
      a_minus_y_full(j1,j2,ell+1) = a - y_u(ell+1);
      y_plus_b_full(j1,j2,ell+1) = y_w(ell+1) + b;
    end
  end
end

% Order zero

Delta = tau2*1i*beta_wp + 1i*beta_up;
M0Inv_11 = tau2*(1i*beta_wp)./Delta;
M0Inv_12 = 1.0./Delta;
M0Inv_21 = -1i*beta_up./Delta;
M0Inv_22 = 1.0./Delta;
Qhat = fft2( eem.*zeta_n(:,:,0+1) );
Rhat = fft2( eem.*psi_n(:,:,0+1) );
uhat = M0Inv_11.*Qhat + M0Inv_12.*Rhat;
what = M0Inv_21.*Qhat + M0Inv_22.*Rhat;

for ell=0:Ny
  un(:,:,ell+1,0+1) = eep.*ifft2( exp(1i*beta_up*y_u(ell+1)).*uhat );
  wn(:,:,ell+1,0+1) = eep.*ifft2( exp(-1i*beta_wp*y_w(ell+1)).*what );
end

% Order n>0

for n=1:N
    
  % Form Qn
  
  Qn = psi_n(:,:,n+1);
  if(n>=1)
    Qn = Qn + ((a-b)/(a*b))*f.*psi_n(:,:,n-1+1);
  end
  if(n>=2)
    Qn = Qn - (1.0/(a*b))*ff.*psi_n(:,:,n-2+1);
  end
  
  % Form Fun, Jun
  
  Fun = zeros(Nx1,Nx2,Ny+1);
  Jun = zeros(Nx1,Nx2);
  
  A1_x1x1 = -(2.0/a)*f_full;
  %A1_x1x2 = 0
  A1_x1y = -(1.0/a)*(a_minus_y_full).*f_x1_full;
  %A1_x2x1 = 0
  A1_x2x2 = -(2.0/a)*f_full;
  A1_x2y = -(1.0/a)*(a_minus_y_full).*f_x2_full;
  A1_yx1 = A1_x1y;
  A1_yx2 = A1_x2y;
  %A1_yy = 0;
  
  A2_x1x1 = (1.0/a^2)*f_full.^2;
  %A2_x1x2 = 0
  A2_x1y = (1.0/a^2)*(a_minus_y_full).*(f_full.*f_x1_full);
  %A2_x2x1 = 0
  A2_x2x2 = (1.0/a^2)*f_full.^2;
  A2_x2y = (1.0/a^2)*(a_minus_y_full).*(f_full.*f_x2_full);
  A2_yx1 = A2_x1y;
  A2_yx2 = A2_x2y;
  A2_yy = (1.0/a^2)*((a_minus_y_full).^2).*(f_x1_full.^2 + f_x2_full.^2);
  
  B1_x1 = -(1.0/a)*f_x1_full;
  B1_x2 = -(1.0/a)*f_x2_full;
  %B1_y = 0;
  
  B2_x1 = (1.0/a^2)*f_full.*f_x1_full;
  B2_x2 = (1.0/a^2)*f_full.*f_x2_full;
  B2_y = (1.0/a^2).*(a_minus_y_full).*(f_x1_full.^2 + f_x2_full.^2);
  
  C1 = -ku2*(2.0/a)*f_full;
  C2 = ku2*(1.0/a^2)*f_full.^2;
  
  if(n>=1)
    u_x1 = dx1p_3d(un(:,:,:,n-1+1),alpha1p,eep,eem,Nx1,Nx2,Ny);
    temp = A1_x1x1.*u_x1;
    Fun = Fun - dx1p_3d(temp,alpha1p,eep,eem,Nx1,Nx2,Ny);
    %A1_x2x1 = 0
    temp = A1_yx1.*u_x1;
    Fun = Fun - dy_3d(temp,Dy,a,Nx1,Nx2,Ny);
    temp = B1_x1.*u_x1;
    Fun = Fun + temp;
    
    u_x2 = dx2p_3d(un(:,:,:,n-1+1),alpha2p,eep,eem,Nx1,Nx2,Ny);
    %A1_x1x2 = 0
    temp = A1_x2x2.*u_x2;
    Fun = Fun - dx2p_3d(temp,alpha2p,eep,eem,Nx1,Nx2,Ny);
    temp = A1_yx2.*u_x2;
    Fun = Fun - dy_3d(temp,Dy,a,Nx1,Nx2,Ny);
    temp = B1_x2.*u_x2;
    Fun = Fun + temp;
    
    u_y = dy_3d(un(:,:,:,n-1+1),Dy,a,Nx1,Nx2,Ny);
    temp = A1_x1y.*u_y;
    Fun = Fun - dx1p_3d(temp,alpha1p,eep,eem,Nx1,Nx2,Ny);
    temp = A1_x2y.*u_y;
    Fun = Fun - dx2p_3d(temp,alpha2p,eep,eem,Nx1,Nx2,Ny);
    %A1_yy = 0
    %B1_y = 0
    
    temp = C1.*un(:,:,:,n-1+1);
    Fun = Fun - temp;
    
    Su = eep.*ifft2( (1i*beta_up).*fft2(eem.*un(:,:,ell_u_top,n-1+1)) );
    Jun = Jun - (1.0/a)*f.*Su;
    
    Qn = Qn - (a/(a*b))*f.*u_y(:,:,ell_u_bottom);
    Qn = Qn + ((a*b)/(a*b))*f_x1.*u_x1(:,:,ell_u_bottom);
    Qn = Qn + ((a*b)/(a*b))*f_x2.*u_x2(:,:,ell_u_bottom);
  end
  
  if(n>=2)
    u_x1 = dx1p_3d(un(:,:,:,n-2+1),alpha1p,eep,eem,Nx1,Nx2,Ny);
    temp = A2_x1x1.*u_x1;
    Fun = Fun - dx1p_3d(temp,alpha1p,eep,eem,Nx1,Nx2,Ny);
    %A2_x2x1 = 0
    temp = A2_yx1.*u_x1;
    Fun = Fun - dy_3d(temp,Dy,a,Nx1,Nx2,Ny);
    temp = B2_x1.*u_x1;
    Fun = Fun + temp;
    
    u_x2 = dx2p_3d(un(:,:,:,n-2+1),alpha2p,eep,eem,Nx1,Nx2,Ny);
    %A2_x1x2 = 0
    temp = A2_x2x2.*u_x2;
    Fun = Fun - dx2p_3d(temp,alpha2p,eep,eem,Nx1,Nx2,Ny);
    temp = A2_yx2.*u_x2;
    Fun = Fun - dy_3d(temp,Dy,a,Nx1,Nx2,Ny);
    temp = B2_x2.*u_x2;
    Fun = Fun + temp;
    
    u_y = dy_3d(un(:,:,:,n-2+1),Dy,a,Nx1,Nx2,Ny);
    temp = A2_x1y.*u_y;
    Fun = Fun - dx1p_3d(temp,alpha1p,eep,eem,Nx1,Nx2,Ny);
    temp = A2_x2y.*u_y;
    Fun = Fun - dx2p_3d(temp,alpha2p,eep,eem,Nx1,Nx2,Ny);
    temp = A2_yy.*u_y;
    Fun = Fun - dy_3d(temp,Dy,a,Nx1,Nx2,Ny);
    temp = B2_y.*u_y;
    Fun = Fun + temp;
    
    temp = C2.*un(:,:,:,n-2+1);
    Fun = Fun - temp;
    
    Qn = Qn + ((a-b)/(a*b))*ff_x1.*u_x1(:,:,ell_u_bottom);
    Qn = Qn + ((a-b)/(a*b))*ff_x2.*u_x2(:,:,ell_u_bottom);
    Qn = Qn - ((a*b)/(a*b))*(f_x1f_x1+f_x2f_x2).*u_y(:,:,ell_u_bottom);
  end
  
  if(n>=3)
    u_x1 = dx1p_3d(un(:,:,:,n-3+1),alpha1p,eep,eem,Nx1,Nx2,Ny);
    Qn = Qn - (1.0/(a*b))*fff_x1.*u_x1(:,:,ell_u_bottom);
    u_x2 = dx2p_3d(un(:,:,:,n-3+1),alpha2p,eep,eem,Nx1,Nx2,Ny);
    Qn = Qn - (1.0/(a*b))*fff_x2.*u_x2(:,:,ell_u_bottom);
    u_y = dy_3d(un(:,:,:,n-3+1),Dy,a,Nx1,Nx2,Ny);
    Qn = Qn - (a/(a*b))*(ff_x1f_x1+ff_x2f_x2).*u_y(:,:,ell_u_bottom);
  end
  
  % Form Fwn, Jwn
  
  Fwn = zeros(Nx1,Nx2,Ny+1);
  Jwn = zeros(Nx1,Nx2);
  
  A1_x1x1 = (2.0/b)*f_full;
  %A1_x1x2 = 0
  A1_x1y = -(1.0/b)*(y_plus_b_full).*f_x1_full;
  %A1_x2x1 = 0
  A1_x2x2 = (2.0/b)*f_full;
  A1_x2y = -(1.0/b)*(y_plus_b_full).*f_x2_full;
  A1_yx1 = A1_x1y;
  A1_yx2 = A1_x2y;
  %A1_yy = 0;
  
  A2_x1x1 = (1.0/b^2)*f_full.^2;
  %A2_x1x2 = 0
  A2_x1y = -(1.0/b^2)*(y_plus_b_full).*(f_full.*f_x1_full);
  %A2_x2x1 = 0
  A2_x2x2 = (1.0/b^2)*f_full.^2;
  A2_x2y = -(1.0/b^2)*(y_plus_b_full).*(f_full.*f_x2_full);
  A2_yx1 = A2_x1y;
  A2_yx2 = A2_x2y;
  A2_yy = (1.0/b^2)*((y_plus_b_full).^2).*(f_x1_full.^2+f_x2_full.^2);
  
  B1_x1 = (1.0/b)*f_x1_full;
  B1_x2 = (1.0/b)*f_x2_full;
  %B1_y = 0;
  
  B2_x1 = (1.0/b^2)*f_full.*f_x1_full;
  B2_x2 = (1.0/b^2)*f_full.*f_x2_full;
  B2_y = -(1.0/b^2).*(y_plus_b_full).*(f_x1_full.^2+f_x2_full.^2);
  
  C1 = kw2*(2.0/b)*f_full;
  C2 = kw2*(1.0/b^2)*f_full.^2;
  
  if(n>=1)
    w_x1 = dx1p_3d(wn(:,:,:,n-1+1),alpha1p,eep,eem,Nx1,Nx2,Ny);
    temp = A1_x1x1.*w_x1;
    Fwn = Fwn - dx1p_3d(temp,alpha1p,eep,eem,Nx1,Nx2,Ny);
    %A1_x2x1 = 0
    temp = A1_yx1.*w_x1;
    Fwn = Fwn - dy_3d(temp,Dy,b,Nx1,Nx2,Ny);
    temp = B1_x1.*w_x1;
    Fwn = Fwn + temp;
    
    w_x2 = dx2p_3d(wn(:,:,:,n-1+1),alpha2p,eep,eem,Nx1,Nx2,Ny);
    %A1_x1x2 = 0
    temp = A1_x2x2.*w_x2;
    Fwn = Fwn - dx2p_3d(temp,alpha2p,eep,eem,Nx1,Nx2,Ny);
    temp = A1_yx2.*w_x2;
    Fwn = Fwn - dy_3d(temp,Dy,b,Nx1,Nx2,Ny);
    temp = B1_x2.*w_x2;
    Fwn = Fwn + temp;
    
    w_y = dy_3d(wn(:,:,:,n-1+1),Dy,b,Nx1,Nx2,Ny);
    temp = A1_x1y.*w_y;
    Fwn = Fwn - dx1p_3d(temp,alpha1p,eep,eem,Nx1,Nx2,Ny);
    temp = A1_x2y.*w_y;
    Fwn = Fwn - dx2p_3d(temp,alpha2p,eep,eem,Nx1,Nx2,Ny);
    %A1_yy = 0
    %B1_y = 0
    
    temp = C1.*wn(:,:,:,n-1+1);
    Fwn = Fwn - temp;
    
    Sw = eep.*ifft2( (-1i*beta_wp).*fft2(eem.*wn(:,:,ell_w_bottom,n-1+1)) );
    Jwn = Jwn + (1.0/b)*f.*Sw;
    
    Qn = Qn - tau2*(b/(a*b))*f.*w_y(:,:,ell_w_top);
    Qn = Qn - tau2*((a*b)/(a*b))*f_x1.*w_x1(:,:,ell_w_top);
    Qn = Qn - tau2*((a*b)/(a*b))*f_x2.*w_x2(:,:,ell_w_top);
  end
  
  if(n>=2)
    w_x1 = dx1p_3d(wn(:,:,:,n-2+1),alpha1p,eep,eem,Nx1,Nx2,Ny);
    temp = A2_x1x1.*w_x1;
    Fwn = Fwn - dx1p_3d(temp,alpha1p,eep,eem,Nx1,Nx2,Ny);
    %A2_x2x1 = 0
    temp = A2_yx1.*w_x1;
    Fwn = Fwn - dy_3d(temp,Dy,b,Nx1,Nx2,Ny);
    temp = B2_x1.*w_x1;
    Fwn = Fwn + temp;
    
    w_x2 = dx2p_3d(wn(:,:,:,n-2+1),alpha2p,eep,eem,Nx1,Nx2,Ny);
    %A2_x1x2 = 0
    temp = A2_x2x2.*w_x2;
    Fwn = Fwn - dx2p_3d(temp,alpha2p,eep,eem,Nx1,Nx2,Ny);
    temp = A2_yx2.*w_x2;
    Fwn = Fwn - dy_3d(temp,Dy,b,Nx1,Nx2,Ny);
    temp = B2_x2.*w_x2;
    Fwn = Fwn + temp;
    
    w_y = dy_3d(wn(:,:,:,n-2+1),Dy,b,Nx1,Nx2,Ny);
    temp = A2_x1y.*w_y;
    Fwn = Fwn - dx1p_3d(temp,alpha1p,eep,eem,Nx1,Nx2,Ny);
    temp = A2_x2y.*w_y;
    Fwn = Fwn - dx2p_3d(temp,alpha2p,eep,eem,Nx1,Nx2,Ny);
    temp = A2_yy.*w_y;
    Fwn = Fwn - dy_3d(temp,Dy,b,Nx1,Nx2,Ny);
    temp = B2_y.*w_y;
    Fwn = Fwn + temp;
    
    temp = C2.*wn(:,:,:,n-2+1);
    Fwn = Fwn - temp;
    
    Qn = Qn - tau2*((a-b)/(a*b))*ff_x1.*w_x1(:,:,ell_w_top);
    Qn = Qn - tau2*((a-b)/(a*b))*ff_x2.*w_x2(:,:,ell_w_top);
    Qn = Qn + tau2*((a*b)/(a*b))*(f_x1f_x1+f_x2f_x2).*w_y(:,:,ell_w_top);
  end
  
  if(n>=3)
    w_x1 = dx1p_3d(wn(:,:,:,n-3+1),alpha1p,eep,eem,Nx1,Nx2,Ny);
    Qn = Qn + tau2*(1.0/(a*b))*fff_x1.*w_x1(:,:,ell_w_top);
    w_x2 = dx2p_3d(wn(:,:,:,n-3+1),alpha2p,eep,eem,Nx1,Nx2,Ny);
    Qn = Qn + tau2*(1.0/(a*b))*fff_x2.*w_x2(:,:,ell_w_top);
    w_y = dy_3d(wn(:,:,:,n-3+1),Dy,b,Nx1,Nx2,Ny);
    Qn = Qn - tau2*(b/(a*b))*(ff_x1f_x1+ff_x2f_x2).*w_y(:,:,ell_w_top);
  end
  
  % Solve elliptic equation
  
  Funhat = zeros(Nx1,Nx2,Ny+1);
  Fwnhat = zeros(Nx1,Nx2,Ny+1);
  for ell=0:Ny
    Funhat(:,:,ell+1) = fft2(eem.*Fun(:,:,ell+1));
    Fwnhat(:,:,ell+1) = fft2(eem.*Fwn(:,:,ell+1));
  end
  Junhat = fft2(eem.*Jun);
  Jwnhat = fft2(eem.*Jwn);
  Qnhat = fft2(eem.*Qn);
  Funhat_p = zeros(Ny+1,1);
  Fwnhat_p = zeros(Ny+1,1);
  
  for j1=1:Nx1
    for j2=1:Nx2
      for ell=0:Ny
        Funhat_p(ell+1) = Funhat(j1,j2,ell+1);
        Fwnhat_p(ell+1) = Fwnhat(j1,j2,ell+1);
      end
      alpha_u = 1.0;
      beta_u = 0.0;
      gamma_u = ku2 - (alpha1p(j1,j2))^2 - (alpha2p(j1,j2))^2;
      alpha_w = 1.0;
      beta_w = 0.0;
      gamma_w = kw2 - (alpha1p(j1,j2))^2 - (alpha2p(j1,j2))^2;
      d_b = -(-1i*beta_wp(j1,j2));
      n_b = 1.0;
      r_b = Jwnhat(j1,j2);
      dw_ab = -1.0;
      du_ab = 1.0;
      rd_ab = zetahat_n(j1,j2,n+1);
      nw_ab = -tau2;
      nu_ab = 1.0;
      rn_ab = Qnhat(j1,j2);
      d_a = -1i*beta_up(j1,j2);
      n_a = 1.0;
      r_a = Junhat(j1,j2);

      [what_p,uhat_p] = solvetwobvp_colloc(Fwnhat_p,alpha_w,beta_w,gamma_w,...
          Funhat_p,alpha_u,beta_u,gamma_u,(2.0/b)*Dy,(2.0/a)*Dy,...
          d_b,n_b,r_b,dw_ab,du_ab,nw_ab,nu_ab,rd_ab,rn_ab,d_a,n_a,r_a);
    
      for ell=0:Ny
        Uhat(j1,j2,ell+1) = uhat_p(ell+1);
        What(j1,j2,ell+1) = what_p(ell+1);
      end
    end
  end
  
  for ell=0:Ny
    un(:,:,ell+1,n+1) = eep.*ifft2(Uhat(:,:,ell+1));
    wn(:,:,ell+1,n+1) = eep.*ifft2(What(:,:,ell+1));
  end

end