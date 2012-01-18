function rbf_svc_demo(n, solver)
%RBF_SVC_DEMO The Demo of RBF kernel SVM for classification
%
%   RBF_SVC_DEMO;
%   RBF_SVC_demo(n);
%
%       Here, n is the number of samples per class. By default it is
%       set to 100.
%
%   RBF_SVC_DEMO(n, solver);
%
%       User can also specifies a customized solver function handle.
%       (default is @mstd_solve using interior-point-convex)
%

%   Created by Dahua Lin, on April 7, 2011
%

%% prepare data

if nargin < 1
    n = 100;
end

if nargin < 2
    solver = svm_default_solver();
end

r0 = 2;
r1 = 12;
r1w = 1;

X0 = randn(2, n) * r0;

t1 = rand(1, n) * (2 * pi);
rho1 = randn(1, n) * r1w + r1;
X1 = [rho1 .* cos(t1); rho1 .* sin(t1)];

%% train svm

X = [X0 X1];
y = [-1 * ones(1, n), ones(1, n)];

C = 10;
sigma = 6;

S = svm_problem('class', X, y, C, {'gauss', sigma});

tic;
R = svm_dual_train(S, [], solver);
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
p = svm_dual_predict(R, gpts);
p = reshape(p, size(xx));
contour(xx, yy, p, [-1, 0, 1]);

% support vectors

hold on;
plot(R.X(1,:), R.X(2,:), 'mo', 'LineWidth', 1, 'MarkerSize', 10);

axis([xm0, xm1, ym0, ym1]);
axis equal;

