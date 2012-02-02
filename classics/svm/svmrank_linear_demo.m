function svmrank_linear_demo(n, solver)
%SVMRANK_LINEAR_DEMO Demo of ranking using linear SVM
%
%   SVMRANK_LINEAR_DEMO(n);
%
%       Here, n is the number of samples. By default it is set to 100.
%
%   SVMRANK_LINEAR_DEMO(n, solver);
%
%       User can also specifies a customized solver function handle.
%       (default is @mstd_solve using interior-point-convex)
%

%   Created by Dahua Lin, on Jan 18, 2012
%       

%% process input

if nargin < 1
    n = 100;
end

if nargin < 2
    solver = svm_default_solver();
end

%% prepare data

A = randn(2, 1);

t = sort(rand(1, n));
t = t(end:-1:1);
X = A * t;

%% train SVM

C = 5;
G = sparse(1:n-1, 2:n, 1, n, n);
S = svm_problem('rank', X, G, C, 'linear');

tic;
w = svm_primal_train(S, solver);
elapsed_t = toc;
fprintf('Training time on %d points = %f sec\n', size(X, 2), elapsed_t);


%% visualize

x = X(1, :);
y = X(2, :);

figure;
scatter(x, y, 400 * (t.^2) );
axis equal;

xlim = get(gca, 'XLim'); x0 = xlim(1); x1 = xlim(2);
ylim = get(gca, 'YLim'); y0 = ylim(1); y1 = ylim(2);

xs = linspace(x0, x1, 200);
ys = linspace(y0, y1, 200);
[xx, yy] = meshgrid(xs, ys);
gpts = [xx(:) yy(:)]';

scores = w' * gpts;
scores = reshape(scores, size(xx));
hold on;
contour(xs, ys, scores, 20);
colorbar;

