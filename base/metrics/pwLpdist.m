% Compute the pairwise Lp-norm distances
%
%   D = pwLpdist(X1, X2, p);
%       computes the Lp-norm distance between pairs of column vectors in 
%       X1 and X2. 
%
%       Suppose the vector dimension is d, then X1 and X2 should be
%       matrices of size d x m and d x n. In this case, the output
%       is a matrix of size m x n, where D(i, j) is the distance between
%       X1(:,i) and X2(:,j).
%

%   Created by Dahua Lin, on Aug 2, 2010
%

