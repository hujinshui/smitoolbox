function mlogireg_demo()
% A simple script to demonstrate the use of mlogireg
%
%   mlogireg_demo;
%

% Created by Dahua Lin, on Jan 25, 2011
%

%% prepare data

n = 300;
m = 3;
d = 6;
a = 2; 
b = 1;

c = randn(2, m) * d;

X = zeros(2, n * m);
for k = 1 : m
    X(:, (k-1) * n + (1:n)) = gen_data(n, c(1,k), c(2,k), a, b, rand() * pi);
end
X = X.';

y = exelem((1:m)', n, 1);

%% regression

f_data = mlogireg(X, y);
f_reg = qreg(1e-3);
f = combfun({f_data, f_reg});

a0 = zeros(2, m);

options = optimset('LargeScale', 'on', 'GradObj', 'on', ...
    'TolX', 1e-8, 'TolFun', 1e-8, 'Display', 'iter');
a = fminunc(f, a0(:), options);
a = reshape(a, 2, m);


%% Visualize data


figure;

% data

colors = {'r', 'g', 'b', 'm'};
for k = 1 : m
    cr = colors{mod(k-1, m) + 1};
    cX = X((k-1) * n + (1:n), :);
    hold on;
    plot(cX(:, 1), cX(:, 2), '.', 'Color', cr);
end

% directions

for k = 1 : m
    hold on;
    plot([0, a(1,k) * 2], [0, a(2,k) * 2], 'Color', 'k', 'LineWidth', 2);
end

% contour

xmin = min(X(:,1));
xmax = max(X(:,1));
ymin = min(X(:,2));
ymax = max(X(:,2));

x0 = xmin - 0.1 * (xmax - xmin);
x1 = xmax + 0.1 * (xmax - xmin);
y0 = ymin - 0.1 * (ymax - ymin);
y1 = ymax + 0.1 * (ymax - ymin);

xx = linspace(x0, x1, 500);
yy = linspace(y0, y1, 500);
[xx, yy] = meshgrid(xx, yy);

Xm = [xx(:), yy(:)];
P = nrmexp(Xm * a, 2);
mp = max(P, [], 2);

hold on;
contour(xx, yy, reshape(mp, size(xx)));


% set range

axis([x0 x1 y0 y1]);
axis equal;

colorbar;


function X = gen_data(n, x0, y0, a, b, t)

X = randn(2, n);
X = bsxfun(@times, [a; b], X);
R = [cos(t) -sin(t); sin(t) cos(t)];
X = bsxfun(@plus, R * X, [x0; y0]);

