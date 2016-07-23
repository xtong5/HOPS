function [u_y] = dy_3d(u,Dy,b,Nx1,Nx2,Ny)

u_y = zeros(Nx1,Nx2,Ny+1);
g = zeros(Ny+1,1);

for j1=1:Nx1
  for j2=1:Nx2
    for ell=0:Ny
      g(ell+1) = u(j1,j2,ell+1);
    end
    u_y(j1,j2,:) = (2.0/b)*Dy*g;
  end
end

return;