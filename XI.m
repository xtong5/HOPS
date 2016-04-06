% hat-function
close all

N=360;
a=0;b=30;c=60;
xi=@(x) (x<b & x >=a).*(x-a)/(b-a)+(x<c & x >=b).*(x-c)/(b-c);
x=2*pi/N*[1:N];
y=xi(x);
polar(x,y)

