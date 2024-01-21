% Define the parameters
N = 100; % Number of grid points
L = 10; % Length of the domain
h = L/N; % Grid spacing
x = linspace(0, L, N+1); % Grid points
y = zeros(N+1, 1); % Initialize the solution vector

% Define the coefficients
a = ones(N+1, 1);
b = -2*ones(N+1, 1);
c = ones(N+1, 1);
d = -h^2*(x.^2+1);

% Apply the boundary conditions
y(1) = 0;
y(N+1) = 0;

% Solve the system of equations
for i = 2:N
    y(i) = (d(i) - a(i)*y(i-1) - c(i)*y(i+1)) / b(i);
end

% Plot the solution
plot(x, y);
xlabel('x');
ylabel('y');
title('Solution of y''''+y''+y+1=0 using finite difference method');
