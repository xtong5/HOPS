function Ln_p=L_n_exterior(u_n_nmo,u_n_nmt,xi_nmo,xi_nmt,Dr_nmo,Dr_nmt,...
    Dp_nmo,Dp_nmt,N_r,N_theta,f,f_theta,a,d,sigma,Y,n)

%compute extra RHS term from BC on interface

u_nmo=zeros(N_theta,1);u_nmt=zeros(N_theta,1);
Dr_u_nmo=zeros(N_theta,1);Dr_u_nmt=zeros(N_theta,1);
Dp_u_nmo=zeros(N_theta,1);Dp_u_nmt=zeros(N_theta,1);

for j = 1:N_theta
        u_nmo(j) = u_n_nmo(j*(N_r+1));
        u_nmt(j) = u_n_nmt(j*(N_r+1));
        Dr_u_nmo(j) = Dr_nmo(j*(N_r+1));
        Dr_u_nmt(j) = Dr_nmt(j*(N_r+1));
        Dp_u_nmo(j) = Dp_nmo(j*(N_r+1));
        Dp_u_nmt(j) = Dp_nmt(j*(N_r+1));
end
if n==1
    Ln = zeros(N_theta,1);
end
if n==2
    Y_u_nmo = ifft(Y.*fft(u_nmo));
    Ln = (d*f.*xi_nmo - a*f.*xi_nmo + sigma*(2*a*d*f.*Dr_u_nmo-d*f_theta.*Dp_u_nmo)...
        +(a*f.*Y_u_nmo-d*f.*Y_u_nmo))/(a*d);
end

if n>2
    Y_u_nmo = ifft(Y.*fft(u_nmo));
    Y_u_nmt = ifft(Y.*fft(u_nmt));
    Ln = (d*f.*xi_nmo - a*f.*xi_nmo -f.*f.*xi_nmt + sigma*(2*a*d*f.*Dr_u_nmo...
        +d*f.*f.*Dr_u_nmt - d*f_theta.*Dp_u_nmo+f.*f_theta.*Dp_u_nmt+d*f_theta.*f_theta.*Dr_u_nmt)...
        +(a*f.*Y_u_nmo-d*f.*Y_u_nmo+f.*f.*Y_u_nmt))/(a*d);
end

Ln_p = fft(Ln);

end