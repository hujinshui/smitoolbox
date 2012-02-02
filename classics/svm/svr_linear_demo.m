function svr_linear_demo(n, solver)
%SVR_LINEAR_DEMO Demo of SVM for linear regression
%
%   SVR_LINEAR_DEMO;
%   SVR_LINEAR_DEMO(n);
%   SVR_LINEAR_DEMO(n, solver);
%
%       Conducts the demo of SVM for linear regression.
%
%       n is the number of samples (default = 200, if omitted).
%       
%       User can also supply a customized solver.
%

% Created by Dahua Lin, on Jan 18, 2011
%

%% process inputs

if nargin < 1
    n = 200;
end

if nargin < 2
    solver = svm_default_solver();
end

%% prepare data

noise_sig = 0.3;

x = rand(1, n) * (2 * pi);

a0 = randn();
b0 = randn();

y = (a0 * x + b0) + randn(1, n) * noise_sig;


%% train SVM

tol = noise_sig;
C = 20;
S = svm_problem('regress', x, y, C, 'linear', tol);
[a, b] = svm_primal_train(S, solver);


%% visualization

figure;
plot(x, y, 'b.');

x0 = 0;
x1 = 2 * pi;
y0 = a * x0 + b;
y1 = a * x1 + b;

hold on;
plot([x0, x1], [y0, y1], 'r-', 'LineWidth', 2);
plot([x0, x1], [y0, y1] + tol, 'm-', 'LineWidth', 1);
plot([x0, x1], [y0, y1] - tol, 'm-', 'LineWidth', 1);

