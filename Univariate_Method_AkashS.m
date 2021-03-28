
%  Non-Linear Optimization
%  Univariate Method by - 
%  AKASH SHAHADE (MT20CDM006)
%  M.Tech., VNIT,Nagpur.

clc
clear
format long

fprintf('Univariate Method by - \n')
fprintf('AKASH SHAHADE - MT20CDM006\n\n');

fprintf('For a function - \n')
fprintf('f = a*X^2 + b*Y^2 + c*X*Y + d*X + e*Y\n\n');
a = input('Enter (a) = ');
b = input('Enter (b) = ');
c = input('Enter (c) = ');
d = input('Enter (d) = ');
e = input('Enter (e) = ');

% Function Definition :
syms X Y;
f = a*X^2 + b*Y^2 + c*X*Y + d*X + e*Y ;

% Initial Guess (Choose Initial Values):
x_cr = input('Enter Initial  X coordinate = ');
y_cr = input('Enter Initial  Y coordinate = ');

x(1) = x_cr;
y(1) = y_cr;
S = [1 0]';
I = [x(1),y(1)]';

% Tolerance and Step-size:
e = 0.01;
i = 1;

% Convergence Parameters:
F = subs(f, [X,Y], [x(1),y(1)]);
f_plus = I + e*S;
F_plus = subs(f, [X,Y], [f_plus(1), f_plus(2)]);
f_minus = I - e*S;
F_minus = subs(f, [X,Y], [f_minus(1), f_minus(2)]);

% Search Direction:
if F_minus < F
    S = -S;
else
    S = S;
end

% Optimization Algorithm:
while ((double(F) > double(F_minus))||(double(F) > double(F_plus)))
    syms h; % Step size
    g = subs(f, [X,Y], [x(i)+S(1)*h,y(i)+h*S(2)]);
    dg_dh = diff(g,h);
    h = solve(dg_dh, h); % Optimal Step Size
    x(i+1) = I(1)+h*S(1); % New x value
    y(i+1) = I(2)+h*S(2); % New y value
    i = i+1;
    I = [x(i),y(i)]'; % Updated Point
    if rem(i,2) == 0
        S = [0 1]';
    else
        S = [1 0]';
    end
    F = subs(f, [X,Y], [x(i),y(i)]);
    f_plus = I + e*S;
    F_plus = subs(f, [X,Y], [f_plus(1), f_plus(2)]);
    f_minus = I - e*S;
    F_minus = subs(f, [X,Y], [f_minus(1), f_minus(2)]);
    if double(F_minus) < double(F)
        S = -S;
    else
        S = S;
    end
end

% Result Table:
Iter = 1:i;
X_coordinate = x';
Y_coordinate = y';
Iterations = Iter';
T = table(Iterations,X_coordinate,Y_coordinate);

% Plots:
fcontour(f,'Fill','On');
hold on;
plot(x,y,'*-r');

% Output:
fprintf('\n Initial Objective Function Value: %d\n\n',subs(f,[X,Y], [x(1),y(1)]));
if ((double(F)>=double(F_minus))||(double(F)>=double(F_plus)))
    fprintf('Minimum succesfully obtained...\n\n');
end
fprintf('Number of Iterations for Convergence: %d\n\n', i);
fprintf('Point of Minima: [%d,%d]\n\n', x(i), y(i));
fprintf('Objective Function Minimum Value after Optimization: %f\n\n', subs(f,[X,Y], [x(i),y(i)]));
disp(T)



