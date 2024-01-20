h=0.03;
k=0.005;
r=k/(h^2);
x=0:h:1;
y=0:k:0.3;
timelevel = length(y);
u = zeros(length(y), length(x));

u(:,1)=0;
u(:, length(x))=0;

for i=1:length(x)
    u(1, i) = sin(pi*x(i));
end

max_iteration = 100;
u = crank_gauss(u, x, timelevel, max_iteration, 1e-6, r);
disp(u);
figure(1);
surf(x, y, u);
xlabel('x');
ylabel('y');
zlabel('z');

figure(2);
[X, Y] = meshgrid(x, y);
ua = exp((-pi^2).*Y).*sin(pi.*X);
surf(x, y, ua); 
 


function u =  crank_gauss(u, x, timelevel, max_iteration, tol, r)
u_old = u;

for k = 1:max_iteration
    for i=2:timelevel
        for j = 2:length(x)-1
              u(i, j)  = (r*u(i-1, j-1) +(2-2*r)*u(i-1, j) +r*u(i-1, j+1) +r*u(i, j-1)+r*u(i, j+1)) / (2+2*r);
        end
    end
    if norm(u-u_old) < tol
        return;
    end
    u_old = u;
  
end
end
  