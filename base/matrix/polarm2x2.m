% Derive the polar representation of a 2 x 2 positive definite matrix
%
%   R = polarm2x2(A);
%       computes the polar representation of a 2x2 positive matrix A.
%       Each polar representation is a column vector with three
%       entries [a; b; theta], such that
%
%           A = R * diag([a b]) * R', 
%
%       where R = [cos(theta), -sin(theta); sin(theta), cos(theta)].
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
%       For any input, the output will be a matrix of size 3 x n, where
%       R(:, i) corresponds to the i-th matrix.
%

%   History
%   -------
%       - Created by Dahua Lin, on June 11, 2010
%

