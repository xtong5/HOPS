function [D_p] = diff_mat_new_interior(p,T_p,Dr,r,k,sigma,a,Z)

%% pre allocation
n = size(r,1);
N_theta = size(p,1);
D_p = zeros(n,n,N_theta);
Diag = diag(r);
A = eye(n); 
DD = Diag*Dr*Diag*Dr + k*k*Diag*Diag;

for i=1:N_theta
    D_p(:,:,i) = DD - p(i)^2*A;
    D_p(1,:,i) = 1/sigma*a*Dr(1,:);
    D_p(1,1,i) = D_p(1,1,i) - Z(i);
    D_p(end,:,i) = Dr(end,:);
    D_p(end,end,i) = D_p(end,end,i)-T_p(i);
    D_p(:,:,i) = D_p(:,:,i)\A;
%     fprintf('Det = %d Cond = %d \n', det(D_p(:,:,i)),cond(D_p(:,:,i)));
end


end