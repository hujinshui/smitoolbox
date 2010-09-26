% Compute pairwise Hamming distances
%
%   D = pwhamdist(X1, X2);
%       compute pairwise Hamming distances between the column vectors in
%       X1 and X2. The hamming distance is defined to be the number of
%       different entries in two vectors.
%
%       Let X1 and X2 be d x m and d x n numerical or logical matrices.
%       Then in the output, D will be a matrix of size m x n, where
%       D(i, j) is the hamming distance between X1(:,i) and X2(:,j).
%       
%       Note that this function supports the following value types:
%       double, single, logical, int32, and uint32.
%

%   Created by Dahua Lin, on Aug 2, 2010
%


