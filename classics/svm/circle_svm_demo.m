function circle_svm_demo(n)
% This function demos to use of kernel SVM on circular data
%
%   circle_svm_demo;
%   circle_svm_demo(n);
%
%       Here, n is the number of samples per class. By default it is
%       set to 100.
%

%   Created by Dahua Lin, on April 7, 2011
%

%% prepare data

if nargin < 1
    n = 100;
end

r0 = 2;
r1 = 9;
r1w = 1;

X0 = randn(2, n) * r0;

t1 = rand(1, n) * (2 * pi);
rho1 = randn(1, n) * r1w + r1;
X1 = [rho1 .* cos(t1); rho1 .* sin(t1)];

%% train svm

X = [X0 X1];
y = [-1 * ones(1, n), ones(1, n)];

C = 10;

kf = @(x1, x2) exp(- pwsqL2dist(x1, x2) / 32 );
K = kf(X, X);
K = (K + K') / 2;
K = adddiag(K, 1e-10);

tic;
svm = kernel_svm.train(X, y, kf, C, 'kermat', K, 'solver', @gurobi_solve);
elapsed_t = toc;
fprintf('Training time on %d samples = %f sec.\n', size(X,2), elapsed_t);

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

% prediction contours

[xx, yy] = meshgrid(linspace(xm0, xm1, 300), linspace(ym0, ym1, 300));
gpts = [xx(:)'; yy(:)'];
p = svm.predict(gpts);
p = reshape(p, size(xx));
contour(xx, yy, p, [-1, 0, 1]);

% support vectors

hold on;
plot(svm.Xs(1,:), svm.Xs(2,:), 'mo', 'LineWidth', 1, 'MarkerSize', 10);

axis([xm0, xm1, ym0, ym1]);
axis equal;

