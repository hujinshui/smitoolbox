function svc_linear_demo(n, solver)
%SVC_LINEAR_DEMO Demo of linear SVM classification
%
%   SVC_LINEAR_DEMO;
%   SVC_LINEAR_DEMO(n);
%   SVC_LINEAR_DEMO(n, solver);
%
%       Here, n is the number of samples per class. By default it is
%       set to 100.
%
%       User can also specifies a customized solver function handle.
%

%   Created by Dahua Lin, on Jan 18, 2012
%

%% process inputs

if nargin < 1
    n = 100;
end

if nargin < 2
    solver = svm_default_solver();
end

%% prepare data

t = rand() * (2 * pi);
tc = t + pi/2 + randn() * 0.5; 
d = 6;

c0 = [cos(tc) sin(tc)] * d;
c1 = -c0;

X0 = gen_data(n, c0(1), c0(2), 3, 1, t);
X1 = gen_data(n, c1(1), c1(2), 3, 1, t);


%% training

X = [X0 X1];
y = [-1 * ones(1, n), ones(1, n)];

C = 5;

S = svm_problem('class', X, y, C);
tic;
[w, b] = svm_primal_train(S, solver);
elapsed_t = toc;
fprintf('Training time on %d points = %f sec\n', size(X, 2), elapsed_t);


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

rgn = [xm0, xm1, ym0, ym1];

% margin lines

impline(w(1), w(2), b-1, rgn, 'Color', 'b');
impline(w(1), w(2), b+1, rgn, 'Color', 'b');

axis([xm0, xm1, ym0, ym1]);
axis equal;


%% Auxiliary functions

function X = gen_data(n, x0, y0, a, b, t)

X = randn(2, n);
X = bsxfun(@times, [a; b], X);
R = [cos(t) -sin(t); sin(t) cos(t)];
X = bsxfun(@plus, R * X, [x0; y0]);


