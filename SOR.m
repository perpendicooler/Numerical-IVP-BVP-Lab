% SOR (Successive Over Relaxation)

n = input('Enter number of equations, n:  ');
A = zeros(n, n + 1);
x1 = zeros(1, n);

A = [9, 1, 3, 7; 1, 12, -4, -15; 3, -4, 20, 10];
x1 = [0, 0, 0];

tol = input('Enter tolerance, tol:  ');
m = input('Enter maximum number of iterations, m:  ');
w = input('Enter the parameter w (omega):  ');

fprintf('\nIteration Table:\n');
fprintf('Iteration\t x1\t\t x2\t\t x3\t\t Error\n');

k = 1;
while k <= m
    err = 0;
    fprintf('%d\t\t', k);
    for i = 1:n
        s = 0;
        for j = 1:n
            s = s - A(i, j) * x1(j);
        end
        s = w * (s + A(i, n + 1)) / A(i, i);
        if abs(s) > err
            err = abs(s);
        end
        x1(i) = x1(i) + s;
        fprintf('%11.8f\t', x1(i));
    end
    
    fprintf('%11.8f\n', err);

    if err <= tol
        break;
    else
        k = k + 1;
    end
end

fprintf('\nThe solution vector after %d iterations is:\n', k);
for i = 1:n
    fprintf('x%d = %11.8f\n', i, x1(i));
end