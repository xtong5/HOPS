function [u_x1] = dx1p_3d(u,alpha1p,eep,eem,Nx1,Nx2,Ny)

u_x1 = zeros(Nx1,Nx2,Ny+1);

for ell=0:Ny
  f = u(:,:,ell+1);
  u_x1(:,:,ell+1) = eep.*ifft2( (1i*alpha1p).*fft2(eem.*f) );
end

return;