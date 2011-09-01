function X = gamma_sample(alpha, beta, n, d)
% Samples from a Gamma distribution
%
%   X = gamma_sample(alpha);       % implicitly set beta = 1
%   X = gamma_sample(alpha, beta);
%   X = gamma_sample(alpha, beta, n);
%   X = gamma_sample(alpha, beta, n, d);
%
%       draws n samples from a gamma distribution with shape parameter 
%       alpha and scale parameter beta.
%
%       Input arguments:
%       - alpha:    can be a scalar or a d x 1 column vector.
%       - beta:     a scalar
%       - n:        the number of samples to draw
%       - d:        the sample space dimension.
%
%       Note that if alpha is a scalar, and d > 1, then d-dimensional
%       samples are drawn, where each component is independent from
%       a Gamma distribution with the same shape param alpha.
%
%       If n is omitted, it is assumed to be 1, and is d is omitted,
%       it is set to size(alpha, 1).
%

% Created by Dahua Lin, on Sep 1, 2011
%

%% verify input

if ~(isfloat(alpha) && ndims(alpha) == 2 && isreal(alpha) && size(alpha,2) == 1)
    error('gamma_sample:invalidarg', ...
        'alpha should be a real scalar or a column vector.');
end
d_ = size(alpha, 1);

if nargin < 2
    beta = 1;
else
    if ~(isfloat(beta) && isscalar(beta) && isreal(beta) && beta > 0)
        error('gamma_sample:invalidarg', ...
            'beta should be a positive real scalar.');
    end
end

if nargin < 3
    n = 1;
else
    if ~(isnumeric(n) && isscalar(n))
        error('gamma_sample:invalidarg', 'n should be a numeric scalar.');
    end
end

if nargin < 4
    d = d_;
else
    if ~(isnumeric(d) && isscalar(d))
        error('gamma_sample:invalidarg', 'd should be a numeric scalar,');
    end
    
    if ~(d_ == 1 || d_ == d)
        error('gamma_sample:invalidarg', 'alpha and d are inconsistent.');
    end
end

%% main

if d == 1
    X = randg(alpha, 1, n);
else
    if d_ == 1
        X = randg(alpha, d, n);
    else
        X = randg(alpha(:, ones(1, n)));
    end
end

if beta ~= 1
    X = X * beta;
end

