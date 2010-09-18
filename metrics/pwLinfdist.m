function D = pwLinfdist(X1, X2)
% Compute the pairwise L_infinity norm distances
%
%   D = pwLinfdist(X1, X2);
%       computes the L_infinity norm distance between pairs of column 
%       vectors in X1 and X2. 
%
%       Suppose the vector dimension is d, then X1 and X2 should be
%       matrices of size d x m and d x n. In this case, the output
%       is a matrix of size m x n, where D(i, j) is the distance between
%       X1(:,i) and X2(:,j).
%

%   Created by Dahua Lin, on Aug 2, 2010
%

%% verify input

if ~(ndims(X1) == 2 && ndims(X2) == 2)
    error('pwL1dist:invalidarg', 'X1 and X2 should be both matrices.');
end


%% main

D = pwLpdist(X1, X2, inf);
    