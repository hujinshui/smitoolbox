function afpclus_demo(n, q)
% Demo of afpclus in clustering 2D points
%
%   afpclus_demo(n);
%   afpclus_demo(n, q);
%
%   Input arguments:
%   - n:        The number of sample points
%   - q:        The quantile ratio (greater q encourages more centers)
%               (default = 0.5)
%

% Created by Dahua Lin, on Mar 27, 2012
%

%% prepare data

if nargin < 2
    q = 0.5;
end

X = rand(2, n);

%% run

D = pwsqL2dist(X);
S = afpsmat(-D, q);

tic;
[M, L] = afpclus(S, 'Display', 'phase');
et = toc;
fprintf('Took %.4f sec.\n\n', et);


%% visualize

figure;

A = sparse(1:n, M(L), true, n, n);
gplot(A, X.');

hold on;
plot(X(1,:), X(2,:), '.');

M = X(:, M);
hold on;
plot(M(1,:), M(2,:), 'ro', 'MarkerSize', 12, 'LineWidth', 2);

