function [Gn] = dno_tfe_helmholtz_lower(un,f,p,alphap,betap,eep,eem,Dy,b,Nx,Ny,N)

Gn = zeros(Nx,N+1);

ell_top = 0 + 1;
f_x = ifft( (1i*p).*fft(f) );

for n=0:N
  u_y = dy(un(:,:,n+1),Dy,b,Nx,Ny);
  Gn(:,n+1) = u_y(:,ell_top);
  if(n>=1)
    u_x = dxp(un(:,:,n-1+1),alphap,eep,eem,Nx,Ny);
    Gn(:,n+1) = Gn(:,n+1) - f_x.*u_x(:,ell_top);
    
    Gn(:,n+1) = Gn(:,n+1) - (1.0/b)*(f.*Gn(:,n-1+1));
  end
  if(n>=2)
    u_x = dxp(un(:,:,n-2+1),alphap,eep,eem,Nx,Ny);
    Gn(:,n+1) = Gn(:,n+1) - (1.0/b)*(f.*(f_x.*u_x(:,ell_top)));

    u_y = dy(un(:,:,n-2+1),Dy,b,Nx,Ny);
    Gn(:,n+1) = Gn(:,n+1) + f_x.*(f_x.*u_y(:,ell_top));
  end
end

return;