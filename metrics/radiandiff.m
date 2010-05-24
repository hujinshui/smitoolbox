function dists = radiandiff(X1, X2)
%RADIANDIFF Computes the radian differences between corresponding vectors
%
%   dists = radiandiff(X1, X2);
%       computes the radian differences between corresponding vectors in X1
%       and X2.
%
%       The radian difference is the angle between two vectors in unit of
%       radian. In mathematics, it can be defined as arccos of normalized
%       correlation, which ranges from 0 (when correlation is 1) to pi
%       (when correlation is -1). If the radian difference between two
%       vectors is pi / 2, they are orthogonal to each other.
%
%       X1 and X2 should be matrices of the same size. Let their size be 
%       d x n, then dists will be a 1 x n vector, with dists(i) being the
%       radian difference betwen X1(:, i) and X2(:, i).
%

%   History
%       - Created by Dahua Lin, on Jun 3, 2008
%

%% parse and verify input arguments

assert(isnumeric(X1) && isnumeric(X2) && ndims(X1) == 2 && ndims(X2) == 2, ...
    'radiandiff:invalidarg', ...
    'X1 and X2 should be both numeric matrices.');

[d, n] = size(X1);
assert(size(X2, 1) == d && size(X2, 2) == n, ...
    'radiandiff:invalidsize', ...
    'X1 and X2 should be of the same size.');


%% main

s1 = sum(X1 .* X1, 1);
s2 = sum(X2 .* X2, 1);

ndp = sum(X1 .* X2, 1) ./ sqrt(s1 .* s2);
ndp(ndp < -1) = -1;
ndp(ndp > 1) = 1;

dists = acos(ndp); 

