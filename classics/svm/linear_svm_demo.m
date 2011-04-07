function linear_svm_demo(n)
% This function demos to use of linear SVM
%
%   linear_svm_demo;
%   linear_svm_demo(n);
%
%       Here, n is the number of samples per class. By default it is
%       set to 50.
%

%   Created by Dahua Lin, on April 7, 2011
%

%% prepare data

if nargin < 1
    n = 50;
end

t = rand() * (2 * pi);
tc = t + pi/2 + randn() * 0.5; 
d = 4;

c0 = [cos(tc) sin(tc)] * d;
c1 = -c0;

X0 = gen_data(n, c0(1), c0(2), 3, 1, t);
X1 = gen_data(n, c1(1), c1(2), 3, 1, t);


%% training

X = [X0 X1];
y = [-1 * ones(1, n), ones(1, n)];

K = X' * X;

C = 10;
P = kernel_svm_prob(K, y, C);
alpha = gurobi_solve(P);

w = bsxfun(@times, X, y) * alpha;
b = - (max(w' * X0) + min(w' * X1)) / 2;



%% visualize

% data points

figure;
plot(X0(1,:), X0(2,:), 'c.'); hold on;
plot(X1(1,:), X1(2,:), 'g.');

% range

xmin = min(X(1,:));
xmax = max(X(1,:));
ymin = min(X(2,:));
ymax = max(X(2,:));

xm0 = xmin - (xmax - xmin) * 0.05;
xm1 = xmax + (xmax - xmin) * 0.05;
ym0 = ymin - (ymax - ymin) * 0.05;
ym1 = ymax + (ymax - ymin) * 0.05;

% margin lines

draw_line(w(1), w(2), b-1, xm0,xm1,ym0,ym1, 'b');
draw_line(w(1), w(2), b+1, xm0,xm1,ym0,ym1, 'b');

% support vectors

sv = X(:, abs(alpha) > 1e-10);
hold on;
plot(sv(1,:), sv(2,:), 'mo', 'LineWidth', 2, 'MarkerSize', 15);


axis([xm0, xm1, ym0, ym1]);
axis equal;


%% Auxiliary functions

function X = gen_data(n, x0, y0, a, b, t)

X = randn(2, n);
X = bsxfun(@times, [a; b], X);
R = [cos(t) -sin(t); sin(t) cos(t)];
X = bsxfun(@plus, R * X, [x0; y0]);

function draw_line(a, b, c, xm0, xm1, ym0, ym1, color)
% draw a line ax + by + c = 0 in range [xm0, xm1, ym0, ym1]

if (abs(a) < abs(b))
    x0 = xm0;
    x1 = xm1;
    y0 = (-a/b) * x0 - (c/b);
    y1 = (-a/b) * x1 - (c/b);
else
    y0 = ym0;
    y1 = ym1;
    x0 = (-b/a) * y0 - (c/a);
    x1 = (-b/a) * y1 - (c/a);
end

line([x0 x1], [y0, y1], 'Color', color);
    


