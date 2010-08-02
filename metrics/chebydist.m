function dists = chebydist(X1, X2)
%CHEBYDIST Computes the Chebyshev distances between corresponding vectors
%
%   dists = chebydist(X1, X2);
%       computes the Chebyshev distances between corresponding vectors in
%       X1 and X2. 
%
%       Chebyshev distance, also known as infinite-norm distance, is
%       defined as
%           
%           d(x, y) = max_i | x(i) - y(i) |
%
%       X1 and X2 should be both d x n matrices, and thus the output dists
%       will be a 1 x n vector, with dists(i) being the distance between
%       X1(:, i) and X2(:, i).
%

%   History
%       - Created by Dahua Lin, on Jun 2, 2008
%

%% parse and verify input arguments

if ~(isnumeric(X1) && isnumeric(X2) && ndims(X1) == 2 && ndims(X2) == 2)
    error('chebydist:invalidarg', ...
        'X1 and X2 should be both numeric matrices.');
end


%% main

dists = max(abs(X1 - X2), [], 1);

