% Finite Difference Method for 1D Heat Conduction
clear;
clc;

% Parameters
L = 1;          % Length of the rod
T = 1;          % Total simulation time
Nx = 100;       % Number of spatial points
Nt = 500;       % Number of time steps
alpha = 0.01;   % Thermal diffusivity

dx = L / (Nx - 1);
dt = T / Nt;

% Initial temperature distribution
u0 = sin(pi * linspace(0, 1, Nx));

% Finite Difference Loop
u = u0;
for n = 1:Nt
    % Interior points (1 to Nx-1)
    for i = 2:Nx-1
        u(i) = u(i) + alpha * dt / dx^2 * (u(i+1) - 2*u(i) + u(i-1));
    end
    
    % Boundary conditions (assuming fixed temperature at both ends)
    u(1) = 0;
    u(Nx) = 0;
    
    % Plot for visualization (optional)

end
    plot(linspace(0, L, Nx), u);
    title(['Time = ', num2str(n*dt)]);
    xlabel('Position');
    ylabel('Temperature');
    ylim([-1, 1]);   % Adjust ylim as needed
    %pause(0.01);     % Pause for visualization
% Final plot
plot(linspace(0, L, Nx), u0, 'b-', linspace(0, L, Nx), u, 'r-');
title('Initial and Final Temperature Distributions');
xlabel('Position');
ylabel('Temperature');
legend('Initial', 'Final');
