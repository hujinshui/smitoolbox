function M = wmedian(X, w)
% Compute the weighted median
%
%   M = wmedian(X, w);
%       computes the weighted median of X (along the first non-singleton
%       dimension of w).
%
%       Let X be an m x n matrix, and w can be either an m x 1 or 1 x n
%       vector. In the former case, it computes the weighted median along 
%       the first dimension, and in the latter case, it computes the
%       weighted median along the second dimension.
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 28, 2010
%

%% verify input

if ~(isfloat(X) && ndims(X) == 2 && isreal(X) && ~issparse(X))
    error('wmedian:invalidarg', 'X should be a non-sparse real matrix.');
end

if ~(isfloat(w) && isvector(w) && isreal(w) && ~issparse(w))
    error('wmedian:invalidarg', 'w should be a non-sparse real vector.');
end

[m, n] = size(X);

if size(w,1) == m && size(w,2) == 1
    dim = 1;        
elseif size(w,1) == 1 && size(w,2) == n
    dim = 2;    
else
    error('wmedian:invalidarg', 'The size of w is inconsistent with X.');
end

if ~isa(w, class(X))
    w = cast(w, class(X));
end

%% main

if dim == 2
    X = X.';
    w = w.';
end

[sX, si] = sort(X, 1);
F = cumsum(w(si), 1);

M = wmedian_cimp(sX, F);

if dim == 2
    M = M.';
end


