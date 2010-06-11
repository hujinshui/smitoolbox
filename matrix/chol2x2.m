% Cholesky decomposition of 2x2 matrices using fast implementation
%
%   R = chol2x2(A);
%       computes the cholesky decomposition of each matrix in A using a 
%       fast implementation. The result is a lower triangle matrix L
%       such that A = L * L'.
%
%       Every matrix in A should be positive definite.
%
%       A can be in either of the following forms:
%       - 2 x 2 matrix
%       - 2 x 2 x n array, with each page being a matrix
%       - 2 x (2 x n) matrix, with all matrices juxtaposed
%       - 4 x n matrix, each column represents a matrix 
%         as [a(1,1), a(2,1), a(1,2), a(2,2)].
%       - 3 x n matrix, each column represents a symmetric matrix
%         as [a(1,1), a(1,2), a(2,2)].
%
%       The output will be in the same form as the input.
%

%   History
%   -------
%       - Created by Dahua Lin, on Apr 7, 2010.
%       - Modified by Dahua Lin, on June 11, 2010.
%
    

