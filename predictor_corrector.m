function [t, y] = predictor_corrector(f, tspan, y0, h)
    % Input:
    %   f: Function handle for the ODE (dy/dt = f(t, y)).
    %   tspan: Time span [t0, tf].
    %   y0: Initial condition.
    %   h: Time step size.

    % Output:
    %   t: Time values.
    %   y: Solution values.
    % Example usage:
%f = @(t, y) -y;  % Example ODE: dy/dt = -y
%tspan = [0, 5];
%y0 = 1;
%h = 0.1;

%[t, y] = predictorCorrectorMethod(f, tspan, y0, h);

% Plot the solution
%plot(t, y);
%xlabel('Time');
%ylabel('Solution y');
%title('Predictor-Corrector Method');


    % Initialize time vector and solution vector
    t = tspan(1):h:tspan(2);
    N = length(t);
    y = zeros(size(t));
    y(1) = y0;

    % Predictor (Adams-Bashforth)
    for i = 2:4
        k1 = h * f(t(i-1), y(i-1));
        k2 = h * f(t(i-1) + h/2, y(i-1) + k1/2);
        k3 = h * f(t(i-1) + h/2, y(i-1) + k2/2);
        k4 = h * f(t(i-1) + h, y(i-1) + k3);
        
        y(i) = y(i-1) + (k1 + 2*k2 + 2*k3 + k4)/6;
    end

    % Corrector (Adams-Moulton)
    for i = 5:N
        % Predictor
        yp = y(i-1) + (h/24)*(55*f(t(i-1),y(i-1)) - 59*f(t(i-2),y(i-2)) + 37*f(t(i-3),y(i-3)) - 9*f(t(i-4),y(i-4)));

        % Corrector
        y(i) = y(i-1) + (h/24)*(9*f(t(i), yp) + 19*f(t(i-1),y(i-1)) - 5*f(t(i-2),y(i-2)) + f(t(i-3),y(i-3)));
    end
end
