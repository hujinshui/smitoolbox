function dists = pwhamdist(X1, X2)
%PWHAMDIST Computes pairwise Hamming distances
%
%   dists = pwhamdist(X1, X2);
%       computes pairwise Hamming distance between the column vectors in X1
%       and X2.
%
%       Hamming distance between two sequence of 0-1 codes is defined as
%
%           d(x, y) = sum_i xor(x(i), y(i))
%
%       which is the number of different code symbols in corresponding
%       positions.
%
%       Suppose X1 and X2 are respectively d x n1 and d x n2 matrices, then
%       dists will be a n1 x n2 matrix, with dists(i, j) being the Hamming
%       distance between X1(:, i) and X2(:, j). Note that both X1 and X2 
%       must be logical matrices.
%
%   dists = pwhamdist(X);
%       computes the pairwise Hamming distances between the column vectors 
%       in X, which is equivalent to pwhamdist(X, X).
%

%   History
%       - Created by Dahua Lin, on Jun 2, 2008
%

%% parse and verify input arguments

if nargin < 2
    X2 = [];
    
    assert(islogical(X1) && ndims(X1) == 2, ...
        'pwhamdist:invalidarg', ...
        'X1 should be a logical matrix.');
    
elseif ~isempty(X2)
    assert(islogical(X1) && islogical(X2) && ndims(X1) == 2 && ndims(X2) == 2, ...
        'pwhamdist:invalidarg', ...
        'X1 and X2 should be both logical matrices.');

    assert(size(X1, 1) == size(X2, 1), ...
        'pwhamdist:invalidarg', ...
        'X1 and X2 should have the same number of rows.');
end

%% main

if isempty(X2)
    X2 = X1;
end

dists = double(pwhamdist_cimp(X1, X2));

