function D = pwL2dist(X1, X2, w)
% Compute the pairwise L2-norm distances
%
%   D = pwL2dist(X1, X2);
%       computes the L2-norm distance between pairs of column vectors in 
%       X1 and X2. 
%
%       Suppose the vector dimension is d, then X1 and X2 should be
%       matrices of size d x m and d x n. In this case, the output
%       is a matrix of size m x n, where D(i, j) is the distance between
%       X1(:,i) and X2(:,j).
%
%   D = pwL2dist(X);
%   D = pwL2dist(X, []);
%       computes the L2-norm distance between pairs of column
%       vectors in X. The implementation for this case is more efficient
%       than pwsqL2dist(X, X), despite that both yield the same result.
%
%   D = pwL2dist(X, [], w);
%   D = pwL2dist(X1, X2, w);
%       computes the weighted L2-norm distance between column vectors
%       in X1 and X2. The weighted L2-norm distance is defined by
%
%           d = sqrt(sum_i w(i) * |x1(i) - x2(i)|^2)
%
%       In the input, w should be a row vector.
%

%   Created by Dahua Lin, on Aug 2, 2010
%


%% main

if nargin < 2
    X2 = [];
end

if nargin < 3
    D = sqrt(pwsqL2dist(X1, X2));
else
    D = sqrt(pwsqL2dist(X2, X2, w));
end

