function lda_binary_demo()
% A simple program to demonstrate the use of lda_binary
%
%   lda_binary_demo;
%

% Created by Dahua Lin, on Nov 22, 2010
%

% prepare model

mu = randn(2, 2) * 4;

sigma = randn(2, 2);
sigma = sigma * sigma';

gd = gaussd.from_mp(mu, gsymat(sigma));

% prepare data

n = 1000;
X = gd.sample(1000);
L = [zeros(1, n), ones(1, n)];

% solve LDA

c0 = sigma \ (mu(:,2) - mu(:,1));
c0 = c0 / norm(c0);  % the truly optimal direction.

c = lda_binary(X, L, 'reg', 1e-8);

fprintf('true optima    = %s\n', num2str(c0.', '%.4f '));
fprintf('solved result  = %s\n', num2str(c.', '%.4f '));

% visualize

figure;
plot(X(1, L==0), X(2, L==0), 'r.');
hold on;
plot(X(1, L==1), X(2, L==1), 'm.');
hold on;
plot(mu(1, :), mu(2, :), 'g+', 'LineWidth', 2, 'MarkerSize', 15);
hold on;
plot_ellipse(gd, 2, 'g-', 'LineWidth', 2);

a0 = mean(mu, 2);
a1 = a0 + c * norm(mu(:,2) - mu(:,1)) / 3;

hold on;
plot([a0(1), a1(1)], [a0(2), a1(2)], 'b-', 'LineWidth', 2);
plot(a0(1), a0(2), 'b.', 'MarkerSize', 20);

axis equal;

