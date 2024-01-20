sx = 0;
sy = 0;
fx = 0.5;
fy = 0.5;
h = 0.125;

k = 0.125;
N = floor((fx - sx)/h)+1;
x = linspace(sx, fx, N);
y = flip(linspace(sy, fy, N));
u = Inf(N, N);
u(N,:) = 0;
u(:,1) = 0;
u(1,:) = 200*x;
u(:,N) = 200*y;
A = zeros((N-2)*(N-2));
B = zeros((N-2)*(N-2), 1);

map = zeros(N);

k = 1;
for i=2:N-1
    for j=2:N-1
        map(i, j) = k;
        k = k+1;
    end
end
k=1;
[A, u, B] = standard_five_point(u, A, B, N, map);
X = zeros((N-2)*(N-2), 1);
fprintf("The generated AX=B matrix using standard five-point method:\n");
disp([array2table(A) array2table(B)]);

X = gaussSeidel(A, B, X, 20, 1e-6);
fprintf("The approximation of internal mesh point using Gauss-Seidel\n");
disp(array2table(X));

X = SOR(A, B, X, 1.25, 20, 1e-6);
fprintf("The approximation of internal mesh point using SOR with w=1.25\n");
disp(array2table(X));

A = zeros((N-2)*(N-2));
B = zeros((N-2)*(N-2), 1);
X = zeros((N-2)*(N-2), 1);

[A, u, B] = five_point_diagonal(u, A, B, N, map);


fprintf("The generated AX=B matrix using diagonal_five_point method:\n");
disp([array2table(A) array2table(B)]);

X = gaussSeidel(A, B, X, 20, 1e-6);
fprintf("The approximation of internal mesh point by diagonal_five_point using Gauss-Seidel\n");
disp(array2table(X));

X = SOR(A, B, X, 1.25, 20, 1e-6);
fprintf("The approximation of internal mesh point by diagonal_five_point  using SOR with w=1.25\n");
disp(array2table(X));


function X = gaussSeidel(A, B, X, maxIterations, tolerance)
    n = length(B);
    iteration = 0;
    error = inf;
    
    while (iteration < maxIterations) && (error > tolerance)
        X_old = X;
        for i = 1:n
            X(i) = (B(i) - A(i,1:i-1)*X(1:i-1) - A(i,i+1:n)*X_old(i+1:n)) / A(i,i);
        end
        iteration = iteration + 1;
        error = norm(X - X_old);
    end
end

function X = SOR(A,B,X,w,maxIterations,tolerance)
    for iteration = 1:maxIterations
        for i = 1:size(A,1)
            X(i) = (1-w)*X(i) + w*(B(i) - (A(i,:)*X - A(i,i)*X(i)))/A(i,i);
        end
        error = norm(A*X - B);
        if error < tolerance
            break;
        end
    end
end

function [A, u, B] = standard_five_point(u, A, B, N, map)
    for i = 2:N-1
        for j = 2:N-1
    %         The following formula of standard five point is used to generate
    %         system of linear equations to solve internal mesh points;
    %          u(i, j) = (u(i-1, j)+u(i+1, j)+u(i, j-1)+u(i, j+1))/4;
            for k = 1: (N-2)*(N-2)
                cr = map(i, j);
                A(cr, map(i, j)) = 4;
                if u(i-1, j) == Inf
                     A(cr, map(i-1, j)) = -1;
                end
                if u(i+1, j) == Inf
                    A(cr , map(i+1, j)) = -1;
                end
                if u(i, j-1) == Inf
                A(cr, map(i, j-1) ) = -1;
                end
                if u(i, j+1) == Inf
                A(cr, map(i, j+1)) = -1;
                end
            end
            cons = 0;
            if u(i-1, j) ~= Inf
                cons = cons + u(i-1, j);
            end
            if u(i+1, j) ~= Inf
                cons = cons + u(i+1, j);
            end
            if u(i, j-1) ~= Inf
                cons = cons + u(i, j-1);
            end
            if u(i, j+1) ~= Inf
                cons = cons + u(i, j+1);
            end
                B(map(i, j)) = cons;
        end
    end
    return;
end

function [A, u, B] = five_point_diagonal(u, A, B, N, map)
    for i = 2:N-1
        for j = 2:N-1
    %         The following formula of five_point_diagonal is used to generate
    %         system of linear equations to solve internal mesh points;
    %       u(i, j) = (u(i-1, j-1) + u(i-1, j+1) + u(i+1, j+1) + u(i+1, j-1) )/4;
            for k = 1: (N-2)*(N-2)
                cr = map(i, j);
                A(cr, map(i, j)) = 4;
                if u(i-1, j-1) == Inf
                     A(cr, map(i-1, j-1)) = -1;
                end
                if u(i-1, j+1) == Inf
                    A(cr , map(i-1, j+1)) = -1;
                end
                if u(i+1, j-1) == Inf
                    A(cr, map(i+1, j-1) ) = -1;
                end
                if u(i+1, j+1) == Inf
                    A(cr, map(i+1, j+1)) = -1;
                end
            end
            cons = 0;
            if u(i-1, j-1) ~= Inf
                cons = cons + u(i-1, j-1);
            end
            if u(i-1, j+1) ~= Inf
                cons = cons + u(i-1, j+1);
            end
            if u(i+1, j-1) ~= Inf
                cons = cons + u(i+1, j-1);
            end
            if u(i+1, j+1) ~= Inf
                cons = cons + u(i+1, j+1);
            end
                B(map(i, j)) = cons;
        end
    end
end