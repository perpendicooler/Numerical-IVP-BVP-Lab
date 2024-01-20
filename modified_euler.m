clear
clc
x0=0; xn=2; h=0.2;
f=@(x,y) (y-x^2+1);
x=x0:h:xn;
y(1)=0.5;
for n= 1:length(x)-1
 yp(n+1) = y(n) + h*f(x(n),y(n));
 y(n+1) = y(n)+ (h/2)*(f(x(n),y(n)) +
f(x(n+1),yp(n+1)));
end
Exact = x.^2+2*x-(exp(x)/2)+1;
absError = abs(Exact-y);
fprintf('\t\txi\t\tExact Value\t\tE''s
Modified\t\tabsError\n')
fprintf('\t\t%.1f\t\t%.8f\t\t%.8f\t\t%.8f\n',[x'
,Exact',y',absError']')
