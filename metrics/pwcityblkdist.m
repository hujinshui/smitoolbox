function dists = pwcityblkdist(X1, X2, w)
%PWCITYBLKDIST Computes pairwise City-block distances
%
%   dists = pwcityblkdist(X1, X2);
%       computes the city-block distances between the columns in X1 and X2.
%
%       The city-block distance, also known as L1 distance, rectilinear 
%       distance, or Manhattan distance, is defined as
%
%           d(x, y) =  sum_i | x(i) - y(i) |
%
%       Suppose X1 and X2 are respectively d x n1 matrix and d x n2 matrix,
%       then dists is a matrix of size n1 x n2, such that dists(i, j) is 
%       the city-block distance between X1(:, i) and X2(:, j).
%
%   dists = pwcityblkdist(X1, X2, w);
%       computes the weighted city-block distances defined as
%
%           d(x, y) = sum_i w(i) | x(i) - y(i) |
%
%       w should be an d x 1 column vector specifying the weights of
%       different components, all weights should be non-negative.
%
%   The following simplified syntax are also supported.
%
%   dists = pwcityblkdist(X);
%   dists = pwcityblkdist(X, []);
%       computes the pairwise city-block distances (or their squares) for
%       the columns in X.
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
        'pwcityblkdist:invalidarg', ...
        'X1 should be a numeric matrix.');
else
    assert(isnumeric(X1) && isnumeric(X2) && ndims(X1) == 2 && ndims(X2) == 2, ...
        'pwcityblkdist:invalidarg', ...
        'X1 and X2 should be both numeric matrices.');
    
    assert(size(X1, 1) == size(X2, 1), ...
        'pwcityblkdist:invalidsize', ...
        'Columns in X1 and X2 differ in length.');
end

if nargin < 3
    weighted = false;
    
else
    d = size(X1, 1);
    assert(isnumeric(w) && ndims(w) == 2 && size(w,1) == d && size(w,2) == 1, ...
        'pwcityblkdist:invalidarg', ...
        'w should be a d x 1 numeric vector.');
    
    assert(all(w >= 0), ...
        'pwcityblkdist:negativeweight', ...
        'all weights should be non-negative.');
    
    weighted = true;
    
end


%% main

if weighted
    X1 = bsxfun(@times, X1, w);
    
    if isempty(X2)
        X2 = X1;
    else
        X2 = bsxfun(@times, X2, w);
    end    
else
    if isempty(X2)
        X2 = X1;
    end
end
        
dists = pwminkdist_cimp(X1, X2, 1);


