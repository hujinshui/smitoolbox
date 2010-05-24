function dists = pwchebydist(X1, X2)
%PWCHEBYDISTS Computes pairwise Chebyshev distances
%
%   dists = pwchebydist(X1, X2);
%       computes pairwise Chebyshev distances between columns in X1 and X2.
%
%       Chebyshev distance, also known as infinite-norm distance, is
%       defined as
%           
%           d(x, y) = max_i | x(i) - y(i) |
%
%       Suppose X1 and X2 are respectively d x n1 and d x n2 matrices, then
%       dists will be an n1 x n2 matrix, with
%
%           dists(i, j) = d(X1(:, i), X2(:, j));
%
%   dists = pwchebydist(X);
%       computes pairwise Chebyshev distances between columns in X, which
%       is equivalent to pwchebydist(X1, X2);
%

%   History
%       - Created by Dahua Lin, on Jun 2, 2008
%

%% parse and verify input arguments

if nargin < 2 
    X2 = [];    
end

if isempty(X2)
    assert(isnumeric(X1) && ndims(X1) == 2, ...
        'pwchebydist:invalidarg', ...
        'X1 should be a numeric matrix.');
else
    assert(isnumeric(X1) && isnumeric(X2) && ndims(X1) == 2 && ndims(X2) == 2, ...
        'pwchebydist:invalidarg', ...
        'X1 and X2 should be both numeric matrices.');
    
    assert(size(X1, 1) == size(X2, 1), ...
        'pwchebydist:invalidsize', ...
        'Columns in X1 and X2 differ in length.');
end

%% main

if isempty(X2)
    X2 = X1;
end

dists = pwminkdist_cimp(X1, X2, inf);