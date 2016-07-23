function [Gn] = dno_tfe_helmholtz_lower_3d(un,f,p1,p2,alpha1p,alpha2p,betap,...
    eep,eem,Dy,b,Nx1,Nx2,Ny,N)

Gn = zeros(Nx1,Nx2,N+1);

ell_top = 0 + 1;
f_x1 = ifft2( (1i*p1).*fft2(f) );
f_x2 = ifft2( (1i*p2).*fft2(f) );

for n=0:N
  u_y = dy_3d(un(:,:,:,n+1),Dy,b,Nx1,Nx2,Ny);
  Gn(:,:,n+1) = u_y(:,:,ell_top);
  if(n>=1)
    u_x1 = dx1p_3d(un(:,:,:,n-1+1),alpha1p,eep,eem,Nx1,Nx2,Ny);
    Gn(:,:,n+1) = Gn(:,:,n+1) - f_x1.*u_x1(:,:,ell_top);
    
    u_x2 = dx2p_3d(un(:,:,:,n-1+1),alpha2p,eep,eem,Nx1,Nx2,Ny);
    Gn(:,:,n+1) = Gn(:,:,n+1) - f_x2.*u_x2(:,:,ell_top);
    
    Gn(:,:,n+1) = Gn(:,:,n+1) - (1.0/b)*(f.*Gn(:,:,n-1+1));
  end
  if(n>=2)
    u_x1 = dx1p_3d(un(:,:,:,n-2+1),alpha1p,eep,eem,Nx1,Nx2,Ny);
    Gn(:,:,n+1) = Gn(:,:,n+1) - (1.0/b)*(f.*(f_x1.*u_x1(:,:,ell_top)));
    
    u_x2 = dx2p_3d(un(:,:,:,n-2+1),alpha2p,eep,eem,Nx1,Nx2,Ny);
    Gn(:,:,n+1) = Gn(:,:,n+1) - (1.0/b)*(f.*(f_x2.*u_x2(:,:,ell_top)));

    u_y = dy_3d(un(:,:,:,n-2+1),Dy,b,Nx1,Nx2,Ny);
    Gn(:,:,n+1) = Gn(:,:,n+1) + f_x1.*(f_x1.*u_y(:,:,ell_top))...
        + f_x2.*(f_x2.*u_y(:,:,ell_top));
  end
end

return;