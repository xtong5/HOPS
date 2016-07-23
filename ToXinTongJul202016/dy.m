function [u_y] = dy(u,Dy,b,Nx,Ny)

u_y = zeros(Nx,Ny+1);
g = zeros(Ny+1,1);

for j=1:Nx
  for ell=0:Ny
    g(ell+1) = u(j,ell+1);
  end
  u_y(j,:) = (2.0/b)*Dy*g;
end

return;