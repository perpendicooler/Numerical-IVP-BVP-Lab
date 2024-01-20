clc
clearvars

f = @(t, y) y - t^2 + 1;

a=0;
b=2;

n = 10;

y(1) = 0.5;
h = (b-a)/(n*1.0);

for i = 1:n;
    ti = a + (i-1)*h;
    k1 = h*f(ti,y(i));
    k2 = h*f(ti + h/2, y(i)+k1/2);
    y(i+1) = y(i) + k2
end

syms t
exact_sol = (t + 1)^2 - 0.5 * exp(t);

for i = 1 : n+1
     ti = a + (i-1)*h;
    exact = subs(exact_sol, ti);
    fprintf("t: %0.2f \ty: %0.7f\t exact: %0.7f \t%0.9f \n", a+(i-1)*h, y(i), exact, abs(y(i)-exact));
end