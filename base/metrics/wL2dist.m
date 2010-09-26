function D = wL2dist(X1, X2, w)
% Compute weighted L2-norm distances between corresponding vectors
%
%   D = wL2dist(X1, X2, w);
%       computes weighted L2-norm distances between the corresponding 
%       column vectors in X1 and X2, defined by
%
%           d = sqrt(sum_i w(i) * |x1(i) - x2(i)|^2)
%
%       X1 and X2 should be matrices of the same size. Suppose their size
%       is d x n, then w should be a matrix of size d x m. Each column of
%       w gives a set of weights. The output is of size m x n, where
%       D(i, :) corresponds to the weights given in w(:, i).
%

%   Created by Dahua Lin, on Aug 2, 2010
%

%% main

D = sqrt(wsqL2dist(X1, X2, w));

