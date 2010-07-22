function dists = mahdist(X1, X2, A)
%MAHDIST Computes the Mahalanobis distances
%
%   dists = mahdist(X1, X2, A);
%       computes the Mahalanobis distances between the corresponding
%       columns in X1 and X2.
%
%       Mahalanobis distances are defined as follows
%
%           d(x, y) = sqrt((x - y)' * A * (x - y));
%
%       In a d-dimensional space, both X1 and X2 should be d x n matrices,
%       and A should be a d x d symmetric matrix.
%
%       In the output, dists will be a 1 x n vector, with 
%           
%           dists(i) = dist(X1(:,i), X2(:,i));
%
%   dists = mahdist(X1, X2, A, 'square');
%       computes the squares of the Mahalanobis distances, i.e. not
%       taking the square root in the above formula.
%
%   Remarks
%       - The matrix A should be a positive definite matrix, the
%         function in itself does not perform the checking for efficiency.
%

%   History
%       - Created by Dahua Lin, on Jun 2, 2008
%       - Modified by Dahua Lin, on Jul 22, 2010
%           - based on sqmahdist.m
%

dists = sqrt(sqmahdist(X1, X2, A));
