function [R, I] = top_k(X, op, K, dim)
% Find the top k elements
%
%   R = top_k(X, 'min', K);
%   R = top_k(X, 'max', K);
%
%       If X is a row/column vector, it returns a row/column vector of
%       length K, which contains the K smallest elements sorted in 
%       ascending order, or the K largest elements sorted in descending
%       order.
%
%       If X is a matrix, it returns a matrix R with K rows, where R(:,i)
%       corresponds to the K top elements for X(:,i).
%
%   R = top_k( ..., dim);
%
%       One can additionally specify the dimension along which the
%       too k elements are extracted.
%
%   [R, I] = top_k( ... );
%       
%       This statement additionally returns the indices of the top-K
%       elements within their own vector.
%
%   Remarks
%   -------
%       1. The implementation is based on a partial sorting algorithm
%          implemented in C++. It is expected to be more efficient 
%          then first launching a full sorting.
%          The complexity is O(n + k logk).
%

%   History
%   -------
%       - Created by Dahua Lin, on Nov 11, 2010
%       - Modified by Dahua Lin, on Mar 28, 2011
%           - re-implemented the core mex using bcslib
%

%% verify input 

if ~(isnumeric(X) && ndims(X) == 2)
    error('top_k:invalidarg', 'X should be a numeric matrix.');
end

if ~(ischar(op))
    error('top_k:invalidarg', 'The 2nd argument should be a string.');
end

switch op
    case 'max'
        code = 1;
    case 'min'
        code = 2;
    otherwise
        error('top_k:invalidarg', 'The 2nd argument is invalid.');
end

if ~(isnumeric(K) && isscalar(K) && K >= 1)
    error('top_k:invalidarg', 'K should be a positive integer scalar.');
end
K = double(K);

if nargin < 4
    dim = 0;
else
    if ~(isnumeric(dim) && isscalar(dim) && (dim == 1 || dim == 2))
        error('top_k:invalidarg', 'dim should be either 1 or 2.');
    end
    dim = double(dim);
end

% check K

if dim == 0  
    if isvector(X)
        check_K(K, numel(X));
    else
        check_K(K, size(X,1));    
    end            
else
    check_K(K, size(X, dim));
end       

%% main

if nargout <= 1
    R = top_k_cimp(X, K, code, dim);
else
    [R, I] = top_k_cimp(X, K, code, dim);
end

            
    
function check_K(K, ub)

if K > ub
    error('top_k:invalidarg', 'The value of K exceeds upper bound.');
end
    
