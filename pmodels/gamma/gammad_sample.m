function X = gammad_sample(alpha, beta, n)
%GAMMAD_SAMPLE Samples from a Gamma distribution
%
%   X = gammad_sample(alpha, beta);
%   X = gammad_sample(alpha, beta, n);
%
%       draws n samples from a gamma distribution with shape parameter 
%       alpha and scale parameter beta.
%
%       Input arguments:
%       - alpha:    can be a scalar or a d x 1 column vector.
%       - beta:     can be a scalar or a d x 1 column vector.
%       - n:        the number of samples to draw
%
%       If n is omitted, it is assumed to be 1, and is d is omitted,
%       it is set to size(alpha, 1).
%
%   X = gammad_sample(alpha, beta, [d, n]);
%       
%       Additionally, specifies the dimension of the samples as d.
%       This syntax is particularly useful when you need to sample
%       from multi-dimensional gamma distributions of which both alpha
%       and beta parameters are scalars.
%

% Created by Dahua Lin, on Sep 1, 2011
% Modified by Dahua Lin, on Dec 26, 2011

%% verify input

if ~(isfloat(alpha) && isreal(alpha) && ndims(alpha) == 2 && size(alpha,2) == 1)
    error('gammad_sample:invalidarg', ...
        'alpha should be a real scalar or a column vector.');
end
if ~(isfloat(beta) && isreal(beta) && ndims(beta) == 2 && size(beta,2) == 1)
    error('gammad_sample:invalidarg', ...
        'beta should be a positive real scalar.');
end

d1 = size(alpha, 1);
d2 = size(beta, 1);

if ~(d1 == 1 || d2 == 1 || d1 == d2)
    error('gammad_sample:invalidarg', ...
        'Inconsistent dimensions between alpha and beta.');
end
d_ = max(d1, d2);

if nargin < 3
    n = 1;
    d = d_;
else
    if ~(isnumeric(n) && (isscalar(n) || numel(n) == 2))
        error('gammad_sample:invalidarg', ...
            'n should be a numeric scalar or pair.');
    end
    
    if isscalar(n)
        d = d_;
    else
        d = n(1);
        n = n(2);
        if ~(d_ == 1 || d_ == d)
            error('gammad_sample:invalidarg', ...
                'The sample dimension is inconsistent.');
        end
    end
end


%% main

if d == 1
    X = randg(alpha, 1, n);
else
    if d1 == 1
        X = randg(alpha, d, n);
    else
        X = randg(alpha(:, ones(1, n)));
    end
end

if isscalar(beta)
    if beta ~= 1
        X = X * beta;
    end
else
    if n == 1
        X = X .* beta;
    else
        X = bsxfun(@times, X, beta);
    end
end



