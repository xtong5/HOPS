function [Qn_tfe] = IIO_new_tfe_interior(Un,Dr_Un,Dp_Un,f,f_theta,a,c,N,sigma,Y)
d=a-c;
Gn_tfe = zeros(size(Un));

for i=1:N+1
    if i == 1
        Gn_tfe(:,1) = a*Dr_Un(:,1);
    end
    if i == 2
        Gn_tfe(:,2) = -f.*(1/a+1/d).*Gn_tfe(:,1)+a*Dr_Un(:,2)+...
            2*f.*Dr_Un(:,1)-f_theta.*Dp_Un(:,1)/a;
    end
    if i > 2
        Gn_tfe(:,i) = -f.*(1/a+1/d).*Gn_tfe(:,i-1)-f.*f.*Gn_tfe(:,i-2)/(a*d)...
        +a*Dr_Un(:,i)+2*f.*Dr_Un(:,i-1)+(f.*f+f_theta.*f_theta).*Dr_Un(:,i-2)/a...
        -f_theta.*Dp_Un(:,i-1)/a-f.*f_theta.*Dp_Un(:,i-2)/(a*d);
    end
end

Qn_tfe = 1/sigma*Gn_tfe - ifft(Y.*fft(Un)); 

end
