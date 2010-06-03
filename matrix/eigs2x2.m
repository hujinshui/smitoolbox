function [e, t] = eigs2x2(A)
% Eigenvalue analysis of 2x2 matrices using fast implementation
%
%   e = eigs2x2(A);
%       computes the eigenvalues of each matrix in A using a 
%       fast implementation. Every matrix in A should be symmetric.
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
%       The output e is a 2 x n matrix, where e(1,:) gives the larger
%       eigenvalue, and e(2,:) gives the smaller eigenvalue. 
%
%   [e, t] = eigs2x2(A);
%       additionally returns the rotation radius t, which can be used
%       to determine the eigenvectors.
%
%       The corresponding eigenvectors are [cos(t) -sin(t)]' and
%       [sin(t) cos(t)]'.
%

% Created by Dahua Lin, on Apr 7, 2010
%

%% verify input arguments

if ~(isfloat(A) && ~issparse(A) && isreal(A))
    error('eigs2x2:invalidarg', ...
        'A should be a non-sparse real-valued array.');
end

[m, n] = size(A);

if ~((m == 2 && mod(n, 2) == 0) || m == 3 || m == 4)
    error('eigs2x2:invalidarg', 'The size of A is invalid.');
end

%% main

if m == 2
    A = reshape(A, 4, []);
    A = A([1 2 4], :);    
elseif m == 4
    A = A([1 2 4], :);
end

R = mat2x2_cimp(A, 6);  

e = R([1 2], :);
t = R(3, :);

