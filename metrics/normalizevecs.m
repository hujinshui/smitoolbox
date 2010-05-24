function Y = normalizevecs(X, p, op)
%NORMALIZEVECS Normalizes the vectors to have unitary norms
%
%   Y = normalizevecs(X, p);
%       computes the normalized vectors with unit p-norm from input vectors.
%
%       Suppose X is a m x n matrix, then Y will be a matrix of the same
%       size. Y(:,i) is proportional to X(:,i) and has unit p-norm.
%
%   Y = normalizevecs(X);
%       computes the normalized vectors with unit L2-norm, which is
%       equivalent to normalizevecs(X, 2);
%
%   Y = normalizevecs(X, p, 'row');
%       computes the normalized row vectors with unit p-norm. In the
%       output, Y(i, :) is proportional to X(i, :) and has unit p-norm.
%

%   History
%       - Created by Dahua Lin, on Jun 3, 2008
%

%% parse and verify input arguments

assert(isnumeric(X) && ndims(X) == 2, 'normalizevecs:invalidarg', ...
    'X should be a numeric matrix.');

if nargin < 2
    p = 2;
    
else
    assert(isscalar(p) && p >= 1, 'normalizevecs:invalidarg', ...
        'p should be a positive scalar with p >= 1.');
end

if nargin < 3
    dim = 1;
else
    assert(strcmp(op, 'row'), 'normalizevecs:invalidarg', ...
        'The 3rd argument can only be ''row'' is specified.');
    
    dim = 2;
end


%% main

if p == 1
    
    Y = bsxfun(@times, X, 1 ./ sum(abs(X), dim));
    
elseif p == 2
    
    Y = bsxfun(@times, X, 1 ./ sqrt(sum(X .* X, dim)) );
    
elseif isinf(p)
    
    Y = bsxfun(@times, X, 1 ./ max(abs(X), [], dim));
    
else
    
    Y = bsxfun(@times, X, 1 ./ (sum(abs(X) .^ p, dim) .^ (1/p)) );
    
end