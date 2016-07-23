function [wtilde,utilde] = solvetwobvp_colloc(ftilde_w,alpha_w,beta_w,gamma_w,...
    ftilde_u,alpha_u,beta_u,gamma_u,D_w,D_u,...
    d_b,n_b,r_b,dw_ab,du_ab,nw_ab,nu_ab,rd_ab,rn_ab,d_a,n_a,r_a)

Ny = length(ftilde_w)-1;
N = Ny+1;

ell_w_min = N;
ell_w_max = 1;
ell_u_min = N;
ell_u_max = 1;

D_w2 = D_w*D_w;
D_u2 = D_u*D_u;

A = zeros(2*N,2*N);
b = zeros(2*N,1);

A(1:N,1:N) = alpha_w*D_w2 + beta_w*D_w + gamma_w*eye(N);
A((N+1):(2*N),(N+1):(2*N)) = alpha_u*D_u2 + beta_u*D_u + gamma_u*eye(N);
b(1:N) = ftilde_w;
b((N+1):(2*N)) = ftilde_u;

% Robin BC at y=-b
ii = ell_w_min;
A(ii,:) = zeros(1,2*N);
A(ii,1:N) = n_b*D_w(ell_w_min,:);
A(ii,ell_w_min) = A(ell_w_min,ell_w_min) + d_b;
b(ii) = r_b;

% Dirichlet BC at y=0
ii = ell_w_max;
A(ii,:) = zeros(1,2*N);
A(ii,ell_w_max) = dw_ab;
A(ii,ell_u_min+N) = du_ab;
b(ii) = rd_ab;

% Neumann BC at y=0
ii = ell_u_min + N;
A(ii,:) = zeros(1,2*N);
A(ii,1:N) = nw_ab*D_w(ell_w_max,:);
A(ii,(N+1):(2*N)) = nu_ab*D_u(ell_u_min,:);
b(ii) = rn_ab;

% Robin BC at y=a
ii = ell_u_max + N;
A(ii,:) = zeros(1,2*N);
A(ii,(N+1):2*N) = n_a*D_u(ell_u_max,:);
A(ii,ell_u_max+N) = A(ell_u_max+N,ell_u_max+N) + d_a;
b(ii) = r_a;

temp = A\b;
wtilde = temp(1:N);
utilde = temp((N+1):2*N);

return;