function dists = minkdist(X1, X2, p, w)
%MINKDIST Computes Minkowski distances
%
%   dists = minkdist(X1, X2, p);
%       computes the Minkowski distances between corresponding columns in
%       X1 and X2.
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
%       For a d-dimensional space, X1 and X2 should be both d x n matrix,
%       then in the output, dists should be a 1 x n vector, with
%
%           dists(i) = d(X1(:, i), X2(:, i));
%
%   dists = minkdist(X1, X2, p, w);
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

assert(isnumeric(X1) && isnumeric(X2) && ndims(X1) == 2 && ndims(X2) == 2, ...
    'minkdist:invalidarg', ...
    'X1 and X2 should be both numeric matrices.');

[d, n] = size(X1);
assert(size(X2, 1) == d && size(X2, 2) == n, ...
    'minkdist:invalidsize', ...
    'X1 and X2 should be of the same size.');

assert(isscalar(p) && isnumeric(p) && p >= 1, ...
    'minkdist:invalidarg', ...
    'p should be a positive scalar with p >= 1.');

if nargin < 4
    weighted = false;
else
    assert(isnumeric(w) && ndims(w) == 2 && size(w,1) == d && size(w,2) == 1, ...
        'minkdist:invalidarg', ...
        'w should be a d x 1 numeric vector.');
    
    assert(all(w >= 0), ...
        'minkdist:negativeweight', ...
        'all weights should be non-negative.');
    
    weighted = true;
end


%% main

if p == 2
    
    D = X1 - X2;
    
    if ~weighted
        dists = sum(D .* D, 1);
    else
        dists = sum(bsxfun(@times, D .* D, w), 1);
    end
    
    dists = sqrt(dists);
    
else
    
    D = abs(X1 - X2);
    
    if p == 1
        if ~weighted
            dists = sum(D, 1);
        else
            dists = sum(bsxfun(@times, D, 3), 1);
        end
        
    elseif isinf(p)
        dists = max(D, [], 1);
        
    else
        if ~weighted
            dists = sum(D .^ p, 1) .^ (1/p);
        else
            dists = sum(bsxfun(@times, D .^ p, w), 1) .^ (1/p);
        end
        
    end
    
end


