function Y = nrmexp(X, dim)
% Compute normalized exponentiation
%
%   Y = nrmexp(X);
%   Y = nrmexp(X, dim);
%       computes normalized exponentiation along the specified dimension.
%       If dim is omitted, the normalization is performed along the 
%       first non-singleton dimension.
%
%       Given a vector x, the normalized exponentiation is a vector y
%       of the same length, which is defined to be
%
%           y_i = exp(x_i) / sum_j exp(x_j).
%
%       To reduce the risk of overflow or underflow, the input values
%       are first appropriately shifted before the computation is
%       performed. This function is useful in converting likelihood 
%       to posterior probabilities.
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 13, 2010
%

%% verify input

if ~isfloat(X)
    error('nrmexp:invalidarg', 'X should be a numeric array.');
end

if nargin < 2
    dim = fnsdim(X);
else
    if ~(isscalar(dim) && dim >= 1 && dim <= ndims(X))
        error('nrmexp:invalidarg', ...
            'dim should be an integer scalar with 1 <= dim <= ndims(X).');
    end
end

%% main

vlen = size(X, dim);

if vlen <= 1
    siz = size(X);
    siz(dim) = 1;
    Y = ones(siz, class(X));
else
    mx = max(X, [], dim);
    Y = exp(bsxfun(@minus, X, mx));
    Y = bsxfun(@times, Y, 1 ./ sum(Y, dim));
end
    
        


