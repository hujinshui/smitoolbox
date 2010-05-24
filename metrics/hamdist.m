function dists = hamdist(X1, X2)
%HAMDIST Computes the hamming distances between corresponding vectors.
%
%   dists = hamdist(X1, X2)
%       computes the hamming distances between corresponding vectors in X1
%       and X2.
%
%       Hamming distance between two sequence of 0-1 codes is defined as
%
%           d(x, y) = sum_i xor(x(i), y(i))
%
%       which is the number of different code symbols in corresponding
%       positions.
%
%       X1 and X2 should be logical matrices with the same size. Let the 
%       size be d x n, then dists will be a 1 x n row vector, with dists(i)
%       being the hamming distance between X1(:,i) and X2(:,i).
%

%   History
%       - Created by Dahua Lin, on Jun 2, 2008
%

%% parse and verify input arguments

assert(islogical(X1) && islogical(X2) && ndims(X1) == 2 && ndims(X2) == 2, ...
    'hamdist:invalidarg', ...
    'X1 and X2 should be both logical matrices.');

[d, n] = size(X1);
assert(size(X2, 1) == d && size(X2, 2) == n, ...
    'hamdist:invalidsize', ...
    'X1 and X2 should be of the same size.');

%% main

dists = sum(X1 ~= X2, 1);
