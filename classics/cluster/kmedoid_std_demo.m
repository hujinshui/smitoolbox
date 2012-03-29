function kmedoid_std_demo(n, K)
% Demo of kmedoid_std in clustering 2D points
%
%   kmedoid_std_demo(n, K);
%
%   Input arguments:
%   - n:        The number of sample points
%   - K:        The number of clusters
%

% Created by Dahua Lin, on Mar 27, 2012
%

%% prepare data

X = rand(2, n);

%% run

cfun = @(i, j) pwsqL2dist(X(:,i), X(:,j));
[M, L] = kmedoid_std({cfun, n}, K, 'Display', 'iter');

assert(isequal(size(M), [1 K]));
assert(isequal(size(L), [1 n]));
assert(all(L >= 1 & L <= K & L == fix(L)));

%% visualize

figure;
plot(X(1,:), X(2,:), '.');
hold on;

M = X(:, M);
plot(M(1,:), M(2,:), 'ro', 'MarkerSize', 12, 'LineWidth', 2);



