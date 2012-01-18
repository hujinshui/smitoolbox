function rbf_svr_demo(n, solver)
%RBF_SVR_DEMO The Demo of RBF kernel SVM for regression
%
%   RBF_SVR_DEMO;
%   RBF_SVR_demo(n);
%
%       Here, n is the number of samples. By default it is set to 200.
%
%   RBF_SVR_DEMO(n, solver);
%
%       User can also specifies a customized solver function handle.
%       (default is @mstd_solve using interior-point-convex)
%

%   Created by Dahua Lin, on April 7, 2011
%


%% process inputs

if nargin < 2
    n = 200;
end

if nargin < 3
    solver = svm_default_solver();
end

%% prepare data

noise_sig = 0.2;

x = sort(rand(1, n) * (2 * pi));
y = sin(x) + randn(1, n) * noise_sig;

%% train SVM

tol = noise_sig;
C = 20;
sigma = 1;
S = svm_problem('regress', x, y, C, {'gauss', sigma}, tol);
R = svm_dual_train(S, [], solver);

% predict

yp = svm_dual_predict(R, x);

%% visualization

figure;
plot(x, y, 'b.');

hold on;
plot(x, yp, 'r-', 'LineWidth', 2);
plot(x, yp + tol, 'm-', 'LineWidth', 1);
plot(x, yp - tol, 'm-', 'LineWidth', 1);


