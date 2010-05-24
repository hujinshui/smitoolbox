function dists = cityblkdist(X1, X2, w)
%CITYBLKDIST Computes the city-block distances between corresponding vectors
%
%   dists = cityblkdist(X1, X2);
%       computes the city-block distances between corresponding vectors in
%       X1 and X2.
%
%       The city-block distance, also known as L1 distance, rectilinear 
%       distance, or Manhattan distance, is defined as
%
%           d(x, y) =  sum_i | x(i) - y(i) |
%
%       X1 and X2 in the input should be of the same size. Suppose they
%       are both d x n matrices, then dists will be a 1 x n vector, with
%
%           dists(i) = d(X1(:, i), X2(:, i));
%
%   dists = cityblkdist(X1, X2, w);
%       computes with weighted city-block distances defined as
%
%           d(x, y) = sum_i w(i) * | x(i) - y(i) |
%
%       w should be an d x 1 column vector specifying the weights of
%       different components, all weights should be non-negative.
%

%   History
%       - Created by Dahua Lin, on Jun 2, 2008
%

%% parse and verify input arguments

assert(isnumeric(X1) && isnumeric(X2) && ndims(X1) == 2 && ndims(X2) == 2, ...
    'cityblkdist:invalidarg', ...
    'X1 and X2 should be both numeric matrices.');

[d, n] = size(X1);
assert(size(X2, 1) == d && size(X2, 2) == n, ...
    'cityblkdist:invalidsize', ...
    'X1 and X2 should be of the same size.');

if nargin < 3
    weighted = false;
else
    assert(isnumeric(w) && ndims(w) == 2 && size(w,1) == d && size(w,2) == 1, ...
        'eucdist:invalidarg', ...
        'w should be a d x 1 numeric vector.');
    
    assert(all(w >= 0), ...
        'eucdist:negativeweight', ...
        'all weights should be non-negative.');
    
    weighted = true;
end


%% main

if ~weighted   
    dists = sum(abs(X1 - X2), 1);
else
    dists = sum(bsxfun(@times, abs(X1 - X2), w), 1);
end

