function D = pwmahdist(X1, X2, A)
%PWMAHDIST Computes the pairwise Mahalanobis distances
%
%   D = pwmahdist(X1, X2, A);
%       computes the pairwise Mahalanobis distances between the columns in
%       X1 and X2.
%
%       Mahalanobis distances are defined as follows
%
%           d(x, y) = sqrt((x - y)' * A * (x - y));
%
%       In a d-dimensional space, both X1 and X2 should be respectively 
%       d x n1 and d x n2 matrices, and A should be a d x d symmetric 
%       positive definite matrix.
%
%       In the output, dists will be a n1 x n2 matrix, with 
%           
%           dists(i, j) = dist(X1(:,i), X2(:,j));
%
%       X2 can be input as an empty matrix, which indicates X2 = X1.
%

%   History
%       - Created by Dahua Lin, on Jun 2, 2008
%       - Modified by Dahua Lin, on Aug 2, 2008
%           - based on pwsqmahdist
%


%% main

D = sqrt(pwsqmahdist(X1, X2, A));


