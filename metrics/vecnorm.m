function v = vecnorm(X, p, dim)
%VECNORM Computes the norm of vectors
%
%   v = vecnorm(X, p);
%       computes the p-norms of the column vectors in X.
%
%       Suppose X is a m x n matrix, then v will be a 1 x n vector, with
%       v(i) being the p-norm of the column vector X(:, i);
%
%       p should be a positive scalar with p >= 1. (p can be infinite)
%
%   v = vecnorm(X);
%       computes the L2-norms of the column vectors in X, which is
%       equivalent to vecnorm(X, 2);
%
%   v = vecnorm(X, p, dim);
%       computes the p-norms of row vectors along dimension dim.
%

%   History
%       - Created by Dahua Lin, on Jun 3, 2008
%

%% parse and verify input arguments

if nargin < 2
    p = 2;
    
else
    if ~(isscalar(p) && p >= 1) 
        error('vecnorm:invalidarg', ...
            'p should be a positive scalar with p >= 1.');
    end
end

if nargin < 3
    dim = 1;
end



%% main

if p == 2
    
    v = sqrt(sum(X .* X, dim));

elseif p == 1
    
    v = sum(abs(X), dim);
    
elseif isinf(p)
    
    v = max(abs(X), [], dim);
    
else
    
    v = sum(abs(X) .^ p, dim) .^ (1/p);
    
end

