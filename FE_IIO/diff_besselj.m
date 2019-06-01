function [diffbesselj]=diff_besselj(p,n,z)
% nth derivative of hankel at 
% wavenumber vector p and scalar z OR
% wavebumber p and vector z
if size(z,1)==1
    diffbesselj=zeros(size(p,1),1);
    for i=0:n
        diffbesselj=diffbesselj+(-1)^i*nchoosek(n,i).*besselj(p.'-n+2*i,z).';   
    end
end
if size(p,1)==1
    diffbesselj=zeros(size(z,1),1);
    for i=0:n
        diffbesselj=diffbesselj+(-1)^i*nchoosek(n,i).*besselj(p-n+2*i,z.').';
    end
end
diffbesselj=diffbesselj./2^n;
end

