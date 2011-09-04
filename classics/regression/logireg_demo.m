function logireg_demo()
% A simple script to demonstrate the use of logireg
%
%   logireg_demo;
%

% Created by Dahua Lin, on Jan 25, 2011
%

%% prepare data

n = 300;
t = rand() * (2 * pi);
tc = t + pi/2 + randn() * 0.5; 
d = 2;

c0 = [cos(tc) sin(tc)] * d;
c1 = -c0;

X0 = gen_data(n, c0(1), c0(2), 3, 1, t);
X1 = gen_data(n, c1(1), c1(2), 3, 1, t);


%% model and solve

X = [X0'; X1'];
y = [zeros(n, 1); ones(n, 1)];

f_data = logireg(X, y);
f_reg = qreg(1e-3);
f = combfun({f_data, f_reg});

options = optimset('LargeScale', 'on', 'GradObj', 'on', ...
    'TolX', 1e-8, 'TolFun', 1e-8, 'Display', 'iter');

a0 = zeros(2, 1);
a = fminunc(f, a0, options);

%% visualize

% data points

figure;
plot(X0(1,:), X0(2,:), 'r.'); hold on;
plot(X1(1,:), X1(2,:), 'g.');

xmin = min(X(:,1));
xmax = max(X(:,1));
ymin = min(X(:,2));
ymax = max(X(:,2));


% direction

hold on;
plot([0, a(1)], [0, a(2)], 'Color', 'b', 'LineWidth', 2);
hold on;
plot(0, 0, '+', 'LineWidth', 2, 'MarkerSize', 10);

% contour

hold on;
xx = linspace(xmin, xmax, 300); 
yy = linspace(ymin, ymax, 300);
[xx, yy] = meshgrid(xx, yy);
pv = 1 ./ (1 + exp(-([xx(:), yy(:)] * a)));

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


