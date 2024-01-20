h = 0.1;
k = 0.1;
tl = 1/k;
xl = 0;
xm = 1;
n = floor((xm - xl)/h)+1;

u=zeros(tl+1, n);
%setting the boundary condition
u(:, 1) = 0.9;
u(:, n) = 0.9;
for i = 1:n
    u(1, i) = 0.9*cos(2*pi*(i-1)*h);
end
alpha = k/h;
ut = zeros(1, n);
for i=1:n
    ut(1, i) = 0;
end

for i=1:tl
    if i==1
        for j = 2:n-1
            u(i+1, j) = (ut(i)*2*k+alpha^2*(u(i, j-1) + u(i, j+1)) + 2*(1-alpha^2)*u(i, j))/2;
        end
        continue;
    end
    for j = 2:n-1
        u(i+1, j) = -u(i-1,j)  + alpha^2*(u(i, j-1) + u(i, j+1)) + 2*(1-alpha^2)*u(i,j);
    end
end