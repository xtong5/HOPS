function [u_x2] = dx2p_3d(u,alpha2p,eep,eem,Nx1,Nx2,Ny)

u_x2 = zeros(Nx1,Nx2,Ny+1);

for ell=0:Ny
  f = u(:,:,ell+1);
  u_x2(:,:,ell+1) = eep.*ifft2( (1i*alpha2p).*fft2(eem.*f) );
end

return;