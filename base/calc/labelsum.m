% Computes the sum respectively for values with different labels
%
%   Y = labelsum(X, K, L);
%       Computes the sum of the values or vectors in X according to the
%       labels given in L.
%
%       In particular, if X is a matrix of size m x n, then L is a vector
%       of integer labels, whose size can be either m x 1 or 1 x n, and
%       the labels in L should be in the range [1, K].
%
%       When L is of size 1 x n, then Y is an m x K matrix, with 
%       Y(:, k) = sum(X(:, L==k), 2); when L is of size m x 1, then 
%       Y is a K x n matrix, with Y(k, :) = sum(X(L == k, :), 1).
%

%   History
%   -------
%       - Created by Dahua Lin, on Jul 9, 2010
%

