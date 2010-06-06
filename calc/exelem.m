function B = exelem(A, nr, nc)
% Construct a new matrix by expanding its elements into sub-matrices
%
%   B = exelem(A, nr, nc);
%       creates a new matrix B by expanding each element in matrix A
%       into a sub-matrix of size nr x nc. 
%
%       Let A is a matrix of size m x n, then in the output, B is a 
%       matrix of the same class as A, whose size is (m x nr) x (n x nc).
%       
%       Example:
%
%           A = [3 4; 5 6];
%           B = repelem(A, 2, 3);
%
%           B = 
%   
%               3 3 3 4 4 4
%               5 5 5 6 6 6
%

%   History
%   -------
%       - Created by Dahua Lin, on Apr 6, 2010
%       - Modified by Dahua Lin, on Jun 6, 2010
%           - change name from repelem to exelem
%           - change error handling (make it lightweight)
%

%% verify input

assert(ndims(A) == 2, 'exelem:invalidarg', ...
    'A should be a matrix with ndims(A) == 2.');

%% main

if isscalar(A)
    B = A(ones(nr, nc));
else
    if nr == 1
        if nc == 1
            B = A;
        else
            J = 1 : size(A, 2);
            J = J(ones(nc, 1), :);
            B = A(:, J(:));
        end
    else
        if nc == 1
            I = 1 : size(A, 1);
            I = I(ones(nr, 1), :);
            B = A(I(:), :);
        else
            I = 1 : size(A, 1);
            I = I(ones(nr, 1), :);
            J = 1 : size(A, 2);
            J = J(ones(nc, 1), :);
            B = A(I(:), :);
            B = B(:, J(:));
        end
    end
end
       