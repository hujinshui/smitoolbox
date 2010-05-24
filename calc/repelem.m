function B = repelem(A, nr, nc)
%REPELEM Construct a new matrix by repeating its elements
%
%   B = repelem(A, nr, nc);
%       creates a new matrix B by repeating the elements in matrix A.
%       In particular, each element in A gives rise to an nr x nc 
%       block in B. If the size of A is m x n, then the size of B is
%       (m x nr) x (n x nc).
%       
%       Example
%       -------
%           A = [3 4; 5 6];
%           B = repelem(A, 2, 3);
%
%           B = 
%   
%               3 3 3 4 4 4
%               5 5 5 6 6 6
%

% Created by Dahua Lin, on Apr 6, 2010
%

%% verify input arguments

assert(ndims(A) == 2, 'repelem:invalidarg', ...
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
       