function B = adddiag(A, a)
% Add values to the diagonal of a matrix
%
%   B = adddiag(A, a);
%       add values in a to the diagonal entries of A, and returns the
%       result.
%
%       Here, A should be a square matrix (say of size n x n), and
%       a can be a scalar or a vector of length n.
%

%% verify input

if ~(isnumeric(A) && ndims(A) == 2)
    error('adddiag:invalidarg', 'A should be a numeric matrix.');
end

[m, n] = size(A);
if m ~= n
    error('adddiag:invalidarg', 'A should be a square matrix.');
end

if size(a, 1) > 1
    a = a.';
end

%% main

if ~issparse(A)  % full matrix
    
    di = 1 + (0:n-1) * (n+1);
    B = A;
    B(di) = B(di) + a;
        
else  % sparse matrix
    
    if isscalar(a)
        B = spdiag(n, a) + A;
    else
        B = spdiag(a) + A;
    end    
    
end
