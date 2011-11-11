function X = dirichlet_sample(K, alpha, n)
% Samples from a Dirichlet distribution
%
%   X = dirichlet_sample(K, alpha, n);
%
%       Draws n samples from a Dirichlet distribution over the 
%       (K-1)-dimensional probability simplex, whose parameter
%       is given by alpha.
%
%       Input arguments:
%       - K:        the sample dimension. Each sample is a vector
%                   of length K that sums to 1.
%
%       - alpha:    the parameter of Dirichlet distribution. It is
%                   either a K x 1 column vector, or a scalar (for
%                   symmetric Dirichlet)
%
%       - n:        the number of samples to draw.
%                   

% Created by Dahua Lin, on Sep 27, 2011
%

%% verify input arguments

if ~(isnumeric(K) && isscalar(K) && K >= 1 && K == fix(K))
    error('dirichlet_sample:invalidarg', ...
        'K should be a positive integer scalar.');
end

if ~( isfloat(alpha) && isreal(alpha) && ...
        (isscalar(alpha) || (ndims(alpha) == 2 && size(alpha,2) ==1)) )
    error('dirichlet_sample:invalidarg', ...
        'alpha should be a real scalar or a real column vector.');
end

if ~(isnumeric(n) && isscalar(n) && n >= 0 && n == fix(n))
    error('dirichlet_sample:invalidarg', ...
        'n should be a non-negative integer scalar.');
end

%% main

if K > 1

    if isscalar(alpha)
        X = randg(alpha, K, n);
    else
        if n == 1
            X = randg(alpha);
        else
            X = randg(alpha(:, ones(1,n)));
        end
    end
    X = bsxfun(@times, X, 1 ./ sum(X, 1));
    
else % K == 1
    X = ones(1, n); 
    
end

