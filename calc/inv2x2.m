function R = inv2x2(A)
% Compute the inverse of 2x2 matrices using fast implementation
%
%   R = inv2x2(A);
%       computes the inverse of each matrix in A using a fast
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
%       The output will be in the same form as the input.
%

% Created by Dahua Lin, on Apr 7, 2010
%

%% verify input arguments

if ~(isfloat(A) && ~issparse(A) && isreal(A))
    error('inv2x2:invalidarg', ...
        'A should be a non-sparse real-valued array.');
end

[m, n] = size(A);

if ~((m == 2 && mod(n, 2) == 0) || m == 3 || m == 4)
    error('inv2x2:invalidarg', 'The size of A is invalid.');
end

%% main

if m == 2 || m == 4
    R = mat2x2_cimp(A, 3);  % for generic matrix
else
    R = mat2x2_cimp(A, 4);  % for compact form of symmetric matrix
end

if size(A, 1) == 2
    R = reshape(R, size(A));
end

