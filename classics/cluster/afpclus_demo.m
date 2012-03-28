function afpclus_demo(n)
% Demo of afpclus in clustering 2D points
%
%   afpclus_demo(n);
%
%   Input arguments:
%   - n:        The number of sample points
%

% Created by Dahua Lin, on Mar 27, 2012
%

%% prepare data

X = rand(2, n);

%% run

S = - pwsqL2dist(X);

ss = S(S < -1e-8);
pref = median(ss);
S(1:(n+1):n^2) = pref;

L = afpclus(S);

M = unique(L);


%% visualize

figure;

A = sparse(1:n, L, true, n, n);
gplot(A, X.');

hold on;
plot(X(1,:), X(2,:), '.');

M = X(:, M);
hold on;
plot(M(1,:), M(2,:), 'ro', 'MarkerSize', 12, 'LineWidth', 2);

