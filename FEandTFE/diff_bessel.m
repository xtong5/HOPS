function [diffbessel]=diff_bessel(c,p,n,z)
% nth derivative of 
% bessel function of first kind c=1
% hankel function of first kind c=2 
% at wavenumber vector p and scalar z OR
% at wavebumber p and vector z
C_nk = ones(n+1,1);
for k = 1:n
    C_nk(k+1) = -C_nk(k)*(n-k+1)/k;
end
if c == 1
% bessel function
   if size(z,1)==1
       diffbessel=zeros(size(p,1),1);
    for i=0:n
        diffbessel=diffbessel+C_nk(i+1).*besselj(p.'-n+2*i,z).';   
    end
   end
   if size(p,1)==1
      diffbessel=zeros(size(z,1),1);
    for i=0:n
        diffbessel=diffbessel+C_nk(i+1).*besselj(p-n+2*i,z').';
    end
   end
else if c == 2
   if size(z,1)==1
       diffbessel=zeros(size(p,1),1);
     for i=0:n
         diffbessel=diffbessel+C_nk(i+1).*besselh(p.'-n+2*i,z).';   
     end
   end
   if size(p,1)==1
       diffbessel=zeros(size(z,1),1);
     for i=0:n
         diffbessel=diffbessel+C_nk(i+1).*besselh(p-n+2*i,z').';
     end
   end
    end
end
diffbessel=diffbessel./2^n;
end


