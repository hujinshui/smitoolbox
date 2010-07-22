function dists = eucdist(X1, X2, w)
% Compute squared Euclidean distances between corresponding vectors
%
%   dists = eucdist(X1, X2);
%       compute the Euclidean distances between the corresponding
%       columns of X1 and X2.
%
%       Suppose both X1 and X2 are matrices of size d x n, then dists
%       will be a 1 x n row vector, with dists(i) being the squared
%       Euclidean distance between X1(:,i) and X2(:,i).
%
%   dists = eucdist(X1, X2, w);
%       compute weighted Euclidean distances between the
%       corresponding vectors of X1 and X2.
%
%       In particular, dists(i) is computed as
%
%           sqrt(\sum_{j=1}^d w(j) * (x1(j) - x2(j))^2)
%
%       Here, x1 and x2 are X1(:,i) and X2(:,i).
%

%   History
%       - Created by Dahua Lin, on Jun 2, 2008
%       - Modified by Dahua Lin, on Jul 22, 2010
%           - based on sqeucdist.m
%


if nargin == 2
    dists = sqrt(sqeucdist(X1, X2));
else
    dists = sqrt(sqeucdist(X1, X2, w));
end



