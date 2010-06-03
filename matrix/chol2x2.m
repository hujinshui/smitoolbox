function R = chol2x2(A)
% Cholesky decomposition of 2x2 matrices using fast implementation
%
%   R = chol2x2(A);
%       computes the cholesky decomposition of each matrix in A using a 
%       fast implementation. The result is a lower triangle matrix.
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

% Created by Dahua Lin, on Apr 7, 2010
%

%% verify input arguments

if ~(isfloat(A) && ~issparse(A) && isreal(A))
    error('chol2x2:invalidarg', ...
        'A should be a non-sparse real-valued array.');
end

[m, n] = size(A);

if ~((m == 2 && mod(n, 2) == 0) || m == 3 || m == 4)
    error('chol2x2:invalidarg', 'The size of A is invalid.');
end

%% main

if m == 2
    siz0 = size(A);
    A = reshape(A, 4, []);
    A = A([1 2 4], :);    
elseif m == 4
    A = A([1 2 4], :);
end
    
R = mat2x2_cimp(A, 5);  

if m ~= 3
    R = [R(1,:); R(2,:); zeros(1, size(R,2)); R(3,:)];
    if m ~= 4
        R = reshape(R, siz0);
    end
end
    

