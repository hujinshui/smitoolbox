function B = rmdiag(A, op)
% Remove the diagonal entries of a square matrix
%
%   B = rmdiag(A);
%   B = rmdiag(A, 'r');
%       
%       Let A be an n x n matrix, then B will be an (n-1) x n matrix, 
%       constructed by removing the diagonal entries, anf lifting
%       the entries which are below the diagonal.
%
%   B = rmdiag(A, 'c');
%
%       With this statement, the output B will be an n x (n-1) matrix,
%       constructed by removing the diagonal entries, and shifting
%       the entries on the right of the diagonal to left.
%
%   Examples
%   --------
%       A = [1 2 3; 4 5 6; 7 8 9]
%
%       rmdiag(A, 'r') => [4 2 3; 7 8 6]
%       rmdiag(A, 'c') => [2 3; 4 6; 7 8]
%

% Created by Dahua Lin, on Nov 13, 2010
%

%% verify input

if ~(ndims(A) == 2 && size(A,1) == size(A,2))
    error('rmdiag:invalidarg', 'A should be a square matrix.');
end

if nargin < 2
    op = 'r';
else
    if ~(isequal(op, 'r') || isequal(op, 'c'))
        error('rmdiag:invalidarg', ...
            'The 2nd argument should be either ''r'' or ''o''.');
    end
end


%% main

n = size(A, 1);

if ~issparse(A)
    
    if op == 'r'
        I = 1 : n*n;
        I(1 + (0:n-1) * (n+1)) = [];
        I = reshape(I, n-1, n);
    else
        I = reshape(reshape(1:n*n, n, n).', 1, n*n);
        I(1 + (0:n-1) * (n+1)) = [];
        I = reshape(I, n-1, n).';
        
    end
    
    B = A(I);
    
else
    [i, j, v] = find(A);
    s0 = find(i == j);
    
    if op == 'r'
        s1 = find(i > j);
        i(s1) = i(s1) - 1;
        m1 = n - 1;
        n1 = n;
    else
        s1 = find(j > i);
        j(s1) = j(s1) - 1;
        m1 = n;
        n1 = n - 1;
    end
    
    i(s0) = [];
    j(s0) = [];
    v(s0) = [];
    
    B = sparse(i, j, v, m1, n1);        
end


