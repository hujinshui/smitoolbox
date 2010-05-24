function v = nrmdot(X1, X2)
%NRMDOT Computes the normalized dot-product between corresponding vectors
%
%   v = nrmdot(X1, X2);
%       computes the normalized dot-product between corresponding column 
%       vectors in X1 and X2.
%
%       Normalized dot-product is defined as the dot-product between
%       vectors after L2-normalization, as
%       
%           normalized dot-product(x, y) = x' * y / (||x|| * ||y||)
%
%       X1 and X2 should be matrices of the same size. Let their size be
%       d x n, then v will be a 1 x n vector with v(i) being the normalized
%       dot product between X1(:, i), and X2(:, i);
%

%   History
%       - Created by Dahua Lin, on Jun 2, 2008
%

%% parse and verify input arguments

assert(isnumeric(X1) && isnumeric(X2) && ndims(X1) == 2 && ndims(X2) == 2, ...
    'nrmdot:invalidarg', ...
    'X1 and X2 should be both numeric matrices.');

[d, n] = size(X1);
assert(size(X2, 1) == d && size(X2, 2) == n, ...
    'nrmdot:invalidsize', ...
    'X1 and X2 should be of the same size.');

%% main

s1 = sum(X1 .* X1, 1);
s2 = sum(X2 .* X2, 1);

v = sum(X1 .* X2, 1) ./ sqrt(s1 .* s2);

v(v > 1) = 1;
v(v < -1) = -1;