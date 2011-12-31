function flda_binary_demo()
% A simple program to demonstrate the use of flda_binary
%
%   flda_binary_demo;
%

% Created by Dahua Lin, on Nov 22, 2010
%

% prepare model

mu = randn(2, 2) * 4;

sigma = randn(2, 2);
sigma = sigma * sigma';

% prepare data

n = 2000;

T = chol(sigma, 'lower');
X1 = bsxfun(@plus, mu(:, 1), T * randn(2, n));
X2 = bsxfun(@plus, mu(:, 2), T * randn(2, n));
X = [X1 X2];
L = [zeros(1, n), ones(1, n)];

% solve LDA

c0 = sigma \ (mu(:,2) - mu(:,1));
c0 = c0 / norm(c0);  % the truly optimal direction.

c = flda_binary(X, L, 'reg', 1e-8);

fprintf('true optima    = %s\n', num2str(c0.', '%.4f '));
fprintf('solved result  = %s\n', num2str(c.', '%.4f '));

% visualize

figure;
plot(X(1, L==0), X(2, L==0), 'g.', 'MarkerSize', 5);
hold on;
plot(X(1, L==1), X(2, L==1), 'b.', 'MarkerSize', 5);
hold on;
plot(mu(1, :), mu(2, :), 'r+', 'LineWidth', 2, 'MarkerSize', 15);

axis equal;
a0 = mean(mu, 2);
impline(c(1), c(2), -c'*a0, [], 'Color', 'r');



