function Normal_n = Nnorm_n(f_theta,N)
% operators by norm of normal vector N

N_theta = size(f_theta,1);
Normal_n = zeros(N_theta,N+1);

Normal_n(:,1) = ones(N_theta,1); %n=0
Normal_n(:,2) = zeros(N_theta,1); %n=1
Normal_n(:,3) = 0.5.*f_theta.*f_theta; %n=2

for n=3:N   
    for m=1:n-1
        Normal_n(:,n+1) = Normal_n(:,n+1)- 0.5.*Normal_n(:,n-m+1).*Normal_n(:,m+1);
    end
end
