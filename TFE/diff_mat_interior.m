function [D_p] = diff_mat_interior(p,T_p,alpha,w,D,x)

%% pre allocation
n = size(x,1);
N_theta = size(p,1);
D_p = zeros(n,n,N_theta);
Diag = diag(x+alpha);
A = eye(n); 
DD = Diag*D*Diag*D + Diag*Diag*w^2;

for i=1:N_theta
    D_p(:,:,i) = DD - p(i)^2*A;
    D_p(1,:,i) = A(1,:);
    D_p(end,:,i)=D(end,:);
    D_p(end,end,i) = D_p(end,end,i)-T_p(i);
%     fprintf('Det = %d Cond = %d \n', det(D_p(:,:,i)),cond(D_p(:,:,i)));
end


end