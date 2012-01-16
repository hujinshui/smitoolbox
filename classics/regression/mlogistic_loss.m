function [v, G] = mlogistic_loss(Z, Y)
%MLOGISTIC_LOSS Multi-class Logistic loss function
%
%   v = MLOGISTIC_LOSS(Z, Y);
%   [v, G] = MLOGISTIC_LOSS(Z, Y);
%
%       Evaluates the multi-class logistic loss function, defined as
%
%           loss(z, y) = - sum_{k=1}^K (y_k * log(p_k))
%
%       with p_k = exp(z_k) / sum_{l=1}^K exp(z_l).
%
%       Here, both z and y are K-dimensional vector. In classification
%       problems, y are typically indicator vectors.
%
%
%       Input arguments:
%       - Z:        The matrix of linear predictors: K x n
%       - Y:        Can be in either of the following forms:
%                   - a K x n matrix of indicator vectors (or soft
%                     assignment vectors). Each column of Y sums to 1.
%                   - a 1 x n vector of class labels. 
%
%       Output arguments:
%       - v:        The loss values [1 x n]
%       - g:        The gradients [K x n]
%

% Created by Dahua Lin, on Jan 1, 2012
%

%% main

[K, n] = size(Z);
zm = max(Z, [], 1);
Z1 = bsxfun(@minus, Z, zm);
E1 = exp(Z1);
se = sum(E1, 1);

if nargout >= 2
    P = bsxfun(@times, E1, 1 ./ se);
end


if size(Y, 1) == 1
    I = Y + (0:n-1) * K;
    v = log(se) - Z1(I);
    if nargout >= 2
        G = P;
        G(I) = G(I) - 1;
    end
    
elseif size(Y, 1) == K
    v = log(se) - sum(Y .* Z1, 1);
    if nargout >= 2
        G = P - Y;
    end
    
else
    error('mlogistic_loss:invalidarg', 'The size of Y is invalid.');
end


