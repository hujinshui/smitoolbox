function dists = sqmahdist(X1, X2, A)
%MAHDIST Computes the squared Mahalanobis distances
%
%   dists = sqmahdist(X1, X2, A);
%       computes the squared Mahalanobis distances between corresponding
%       columns in X1 and X2.
%
%       Squared Mahalanobis distances are defined as follows
%
%           d(x, y) = (x - y)' * A * (x - y);
%
%       In a d-dimensional space, both X1 and X2 should be d x n matrices,
%       and A should be a d x d symmetric matrix.
%
%       In the output, dists will be a 1 x n vector, with 
%           
%           dists(i) = d(X1(:,i), X2(:,i));
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

%% verify input arguments

if ~(ndims(X1) == 2 && isreal(X1) && ndims(X2) == 2 && isreal(X2))
    error('sqmahdist:invalidarg', ...
        'X1 and X2 should be real matrices.');
end

%% main

D = X1 - X2;
dists = sum(D .* (A * D), 1);
