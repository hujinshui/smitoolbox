function kmeans_std_demo(n, K)
% Demo of kmeans_std in clustering 2D points
%
%   kmeans_std_demo(n, K);
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

[M, L] = kmeans_std(X, K, 'Display', 'iter');

assert(isequal(size(M), [2 K]));
assert(isequal(size(L), [1 n]));
assert(all(L >= 1 & L <= K & L == fix(L)));

%% visualize

figure;
plot(X(1,:), X(2,:), '.');
hold on;
plot(M(1,:), M(2,:), 'r+', 'MarkerSize', 12, 'LineWidth', 2);



