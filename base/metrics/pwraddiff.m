function dists = pwraddiff(X1, X2)
%Compute the pairwise radian differences
%
%   dists = pwraddiff(X);
%   dists = pwraddiff(X1, X2);
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
    X2 = [];
end

dists = acos(pwnrmdot(X1, X2));

