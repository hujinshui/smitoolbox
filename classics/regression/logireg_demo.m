function logireg_demo()
% A simple script to demonstrate the use of logireg
%
%   logireg_demo;
%

% Created by Dahua Lin, on Jan 25, 2011
%

%% prepare data

n = 1000;
t = rand() * (2 * pi);
tc = t + pi/2 + randn() * 0.5; 
d = 3;

xc = randn() * 10;
yc = randn() * 10;
dir = [cos(tc) sin(tc)] * d;

X0 = gen_data(n, xc + dir(1), yc + dir(2), 3, 1, t);
X1 = gen_data(n, xc - dir(1), yc - dir(2), 3, 1, t);


%% model and solve

X = [X0, X1];
y = [zeros(1, n), ones(1, n)];

f = logiregf(X, y, [], 1e-3);

a0 = zeros(2+1, 1);
a = bfgsfmin(f, a0, ...
    'MaxIter', 300, 'TolFun', 1e-8, 'TolX', 1e-8, 'Display', 'iter');

fprintf('Boundary: %.4f x + %.4f y + %0.4f = 0\n', a(1), a(2), a(3));


%% visualize

% data points

figure;
plot(X0(1,:), X0(2,:), 'b.', 'MarkerSize', 5); hold on;
plot(X1(1,:), X1(2,:), 'r.', 'MarkerSize', 5);

xmin = min(X(1,:));
xmax = max(X(1,:));
ymin = min(X(2,:));
ymax = max(X(2,:));

% contour

hold on;
xx = linspace(xmin, xmax, 300); 
yy = linspace(ymin, ymax, 300);
[xx, yy] = meshgrid(xx, yy);
zz = [xx(:), yy(:)];
pv = 1 ./ (1 + exp(-(zz * a(1:2) + a(3))));

hold on;
contour(xx, yy, reshape(pv, size(xx)), 0.1:0.1:0.9);
colorbar;

axis([xmin, xmax, ymin, ymax]);
axis equal;




function X = gen_data(n, x0, y0, a, b, t)

X = randn(2, n);
X = bsxfun(@times, [a; b], X);
R = [cos(t) -sin(t); sin(t) cos(t)];
X = bsxfun(@plus, R * X, [x0; y0]);


