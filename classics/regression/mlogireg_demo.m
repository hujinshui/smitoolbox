function mlogireg_demo()
%MLOGIREG_DEMO Demo of Multi-class logistic regression
%
%   mlogireg_demo;
%   
%       This is a simple program that demonstrates the use of 
%       multi-class logistic regression (mlogiregf)
%

%   History
%   -------
%       - Created by Dahua Lin, on Jan 25, 2011
%       - Modified by Dahua Lin, on Jan 15, 2012
%

%% prepare data

n = 300;
K = 3;
d = 6;
a = 2; 
b = 1;

c = randn(2, K) * d;

X = zeros(2, n * K);
for k = 1 : K
    X(:, (k-1) * n + (1:n)) = gen_data(n, c(1,k), c(2,k), a, b, rand() * pi);
end

y = exelem((1:K), 1, n);

%% regression

f = mlogiregf(X, y, [], 1e-3);

a0 = zeros(3, K);
a = bfgsfmin(f, a0(:), ...
    'MaxIter', 300, 'TolFun', 1e-8, 'TolX', 1e-8, 'Display', 'iter');
a = reshape(a, 3, K);


%% Visualize data


figure;

% data

colors = {'r', 'g', 'b', 'm'};
for k = 1 : K
    cr = colors{mod(k-1, K) + 1};
    cX = X(:, (k-1) * n + (1:n));
    hold on;
    plot(cX(1, :), cX(2, :), '.', 'Color', cr);
end

% contour

xmin = min(X(1,:));
xmax = max(X(1,:));
ymin = min(X(2,:));
ymax = max(X(2,:));

x0 = xmin - 0.1 * (xmax - xmin);
x1 = xmax + 0.1 * (xmax - xmin);
y0 = ymin - 0.1 * (ymax - ymin);
y1 = ymax + 0.1 * (ymax - ymin);

xx = linspace(x0, x1, 500);
yy = linspace(y0, y1, 500);
[xx, yy] = meshgrid(xx, yy);
nn = numel(xx);

Xm = [xx(:), yy(:), ones(nn, 1)]';
P = nrmexp(a' * Xm, 1);
mp = max(P, [], 1);

hold on;
contour(xx, yy, reshape(mp, size(xx)), 0.3:0.1:0.9);


% set range

axis([x0 x1 y0 y1]);
axis equal;

colorbar;


function X = gen_data(n, x0, y0, a, b, t)

X = randn(2, n);
X = bsxfun(@times, [a; b], X);
R = [cos(t) -sin(t); sin(t) cos(t)];
X = bsxfun(@plus, R * X, [x0; y0]);

