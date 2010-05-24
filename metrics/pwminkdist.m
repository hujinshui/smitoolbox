function dists = pwminkdist(X1, X2, p, w)
%PWMINKDIST Computes pairwise Minkowski distances
%
%   dists = pwminkdist(X1, X2, p); 
%       computes pairwise Minkowski distances between the columns in X1 and
%       X2.
%
%       Minkowski-distance, also known as p-norm distance, is defined as
%
%           d(x, y) = ( sum_i ( x(i) - y(i) )^p )^(1/p)
%
%       The distance is controlled by a parameter p, which should be a 
%       positive scalar with p >= 1.
%
%       Typical p values include 1, 2, inf, which respectively correspond
%       to city-block distance, Euclidean distance, and Chebyshev distance.
%
%       Suppose X1 and X2 are respectively d x n1 matrix and d x n2 matrix,
%       then dists is a matrix of size n1 x n2, such that dists(i, j) is 
%       the minkowski distance between X1(:, i) and X2(:, j).
%
%   dists = pwminkdist(X1, X2, p, w);
%       computes the weighted Minkowski distance between corresponding 
%       columns in X1 and X2.
%
%       The weighted Minkowski distance is a simple extension, as
%
%           d(x, y) = ( sum_i w(i) * ( x(i) - y(i) )^p )^(1/p)
%
%       The weights are given in the argument w, which should be a
%       d x 1 vector.
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
        'pwminkdist:invalidarg', ...
        'X1 should be a numeric matrix.');
else
    assert(isnumeric(X1) && isnumeric(X2) && ndims(X1) == 2 && ndims(X2) == 2, ...
        'pwminkdist:invalidarg', ...
        'X1 and X2 should be both numeric matrices.');
    
    assert(size(X1, 1) == size(X2, 1), ...
        'pwminkdist:invalidsize', ...
        'Columns in X1 and X2 differ in length.');
end

assert(isscalar(p) && isnumeric(p) && p >= 1, ...
    'pwminkdist:invalidarg', ...
    'p should be a positive scalar with p >= 1.');

if nargin < 4
    weighted = false;
else
    d = size(X1, 1);
    assert(isnumeric(w) && ndims(w) == 2 && size(w,1) == d && size(w,2) == 1, ...
        'pwminkdist:invalidarg', ...
        'w should be a d x 1 numeric vector.');
    
    assert(all(w >= 0), ...
        'pwminkdist:negativeweight', ...
        'all weights should be non-negative.');
    
    weighted = true;
end

%% main

if weighted && ~isinf(p)
    if p == 1
        sw = w;
    elseif p == 2
        sw = sqrt(w);
    else
        sw = w .^ (1 / p);
    end
    
    X1 = bsxfun(@times, X1, sw);
    
    if isempty(X2)
        X2 = X1;
    else
        X2 = bsxfun(@times, X2, sw);
    end    
else
    if isempty(X2)
        X2 = X1;
    end
end


if p == 2
    dists = pweucdist(X1, X2);
    
else
    dists = pwminkdist_cimp(X1, X2, double(p));
    
end


