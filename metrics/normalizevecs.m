function Y = normalizevecs(X, p, dim)
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
%   Y = normalizevecs(X, p, dim);
%       normalizes the vectors along the dimension specified by dim.
%       When dim is omitted, it is equivalent to setting dim to 1.
%

%   History
%       - Created by Dahua Lin, on Jun 3, 2008
%       - Modified by Dahua Lin, on Jul 22, 2010
%           - support user-specified along dimension
%           - simplify error handling
%

%% parse and verify input arguments

if nargin < 2
    p = 2;    
else
    if ~(isscalar(p) && p >= 1) 
        error('normalizevecs:invalidarg', ...
            'p should be a positive scalar with p >= 1.');
    end
end

if nargin < 3
    dim = 1;
end


%% main

Y = bsxfun(@times, X, 1 ./ vecnorm(X, p, dim));
