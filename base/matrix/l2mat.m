function A = l2mat(K, L, w, op)
%L2MAT Assignment matrix from labels
%
%   A = L2MAT(K, L, w);
%       
%       Constructs an assignment matrix A based on a label vector L, such 
%       that 
%
%       If L is a row vector, then A is a K x n matrix, with 
%
%           A(L(i), i) = w(i).
%
%       If L is a column vector, then A is an n x K matrix, with
%
%           A(i, L(i)) = w(i).
%
%
%       Input arguments:
%       - K:        The number of distinct labels. 
%
%       - L:        The vector of labels. The value of L(i) is expected 
%                   in {1, ..., K}. If L(i) is out of this range, it will 
%                   be ignored.       
%
%       - w:        The weights, which can be either a scalar (all weights
%                   being the same) of a vector of length n.
%
%                   The class of A will be the same as that of w.
%
%   A = L2MAT(K, L, w, 'sparse');
%
%       Constructs A as a sparse matrix.
%

%   History
%   -------
%       - Created by Dahua Lin, on Oct 31, 2009
%       - Modified by Dahua Lin, on Apr 7, 2010
%           - change name to l2mat
%           - add support of options: logical and sparse.
%       - Modified by Dahua Lin, on Feb 3, 2010
%           - support user-supplied weights.
%

%% verify input arguments

if ~(isscalar(K) && isnumeric(K) && K == fix(K) && K >= 1)
    error('l2mat:invalidarg', 'K should be a positive integer.');
end

if ~(isnumeric(L) && isvector(L))
    error('l2mat:invalidarg', 'L should be a numeric vector.');
end
n = numel(L);

if ~( (isnumeric(w) || islogical(w)) && ...
        (isscalar(w) || (isvector(w) && numel(w) == n)) )
    error('l2mat:invalidarg', 'w should be a vector of length n.');
end

if nargin >= 4
    if ~strcmp(op, 'sparse')
        error('l2mat:invalidarg', 'The 4th argument is invalid.');
    end
    use_sparse = 1;
else
    use_sparse = 0;
end

%% main

if size(L, 2) == 1
    cf = 1;
else
    cf = 0;
    L = L.';
end

valid = (L >= 1 & L <= K);
if all(valid)
    I = (1:n).';
    J = double(L);
else
    I = find(valid);
    J = L(valid);
    if ~isscalar(w)
        w = w(valid);
    end
end

if size(w, 2) > 1
    w = w.';
end


if use_sparse            
    if cf
        A = sparse(I, J, w, n, K);
    else
        A = sparse(J, I, w, K, n);
    end

else
    if cf
        if islogical(w)
            A = false(n, K);
        else
            A = zeros(n, K, class(w));
        end
        A(I + (J - 1) * n) = w;
    else
        if islogical(w)
            A = false(K, n);
        else
            A = zeros(K, n, class(w));
        end
        A(J + (I - 1) * K) = w;
    end

end


