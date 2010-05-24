function dists = pwradiandiff(X1, X2)
%PWRADIANDIFF Computes the pairwise radian differences
%
%   dists = pwradiandiff(X1, X2);
%       computes the pairwise radian differences between column vectors in
%       X1 and X2.
%
%       The radian difference is the angle between two vectors in unit of
%       radian. In mathematics, it can be defined as arccos of normalized
%       correlation, which ranges from 0 (when correlation is 1) to pi
%       (when correlation is -1). If the radian difference between two
%       vectors is pi / 2, they are orthogonal to each other.
%
%       Suppose X1 and X2 are respectively d x n1 and d x n2 matrices, 
%       then dists will be an n1 x n2 matrix with dists(i, j) being the
%       radian difference between the two vectors X1(:, i) and X2(:, j).
%

%   History
%       - Created by Dahua Lin, on Jun 3, 2008
%

%% main

if nargin < 2
    assert(ndims(X1) == 2, 'pwradiandiff:invalidarg', ...
        'X1 should be a matrix.');
    
    X2 = [];
    
elseif ~isempty(X2)    
    assert(ndims(X1) == 2 && ndims(X2) == 2 && size(X1,1) == size(X2,1), ...
        'pwradiandiff:invalidarg', ...
        'X1 and X2 should be both matrices with the same number of rows.');
    
end

%% main

dists = acos(nrmdot(X1, X2));

