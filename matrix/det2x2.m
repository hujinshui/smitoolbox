% Compute the determinant of 2x2 matrices using fast implementation
%
%   R = det2x2(A);
%       computes the determinant of each matrix in A using a fast
%       implementation.
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
%       The output will be a 1 x n row vector.
%

%   History
%   -------
%       - Created by Dahua Lin, on Apr 7, 2010
%       - Modified by Dahua Lin, on June 20, 2010
%           - based on pure C++ mex implementation
%

