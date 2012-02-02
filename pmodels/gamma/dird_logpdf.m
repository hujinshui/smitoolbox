function L = dird_logpdf(alpha, X, c)
%DIRD_LOGPDF Evaluates log-pdf of Dirichlet distribution
%
%   L = dird_logpdf(alpha, X)
%   L = dird_logpdf(alpha, X, 0);
%   L = dird_logpdf(alpha, X, c);
%
%       Evaluates the log probability density values at given samples
%       with respect to given Dirichlet distribution(s).
%
%       Input arguments:
%       - alpha:        the parameter(s) of the Dirichlet distribution
%                       The size of alpha can be either 1 x 1, d x 1,
%                       1 x m, or d x m.
%
%                       Here, d is the space dimension, and m is the 
%                       number of distributions. If m > 1, then the log
%                       pdf values w.r.t. all distributions will be 
%                       evaluated.
%
%       - X:            the sample matrix, of size d x n. Each column of
%                       X is a sample, which sums to 1.
%
%       - c:            The constant given by multivariate beta function,
%                       which can be pre-computed as mvbetaln(alpha, d).
%
%                       If c is not supplied, the function invokes 
%                       mvbetaln to evaluate it.
%
%                       If c is set to zero, then only the linear terms
%                       are computed, without adding the constant.
%
%       Output arguments:
%       - L:            the log pdf values. The size of L is m x n, where
%                       L(k, i) is the log-pdf at X(:,i) with respect to 
%                       the k-th parameter in alpha.
%

% Created by Dahua Lin, on Dec 26, 2011
%

%% verify input arguments

if ~(isfloat(alpha) && isreal(alpha) && ndims(alpha) == 2)
    error('dird_logpdf:invalidarg', 'alpha should be a real matrix.');
end
[da, m] = size(alpha);

if ~(isfloat(X) && isreal(X) && ndims(X) == 2)
    error('dird_logpdf:invalidarg', 'X should be a real matrix.');
end
d = size(X, 1);

if ~(da == 1 || da == d)
    error('dird_logpdf:invalidarg', 'alpha and X have inconsistent dimensions.');
end

if nargin >= 3
    if ~( isequal(c, 0) || ...
            (isfloat(c) && isreal(c) && isequal(size(c), [1 m])) )
        error('gammad_logpdf:invalidarg', ...
            'c should be either zero or a 1 x m real vector.');
    end
    calc_c = 0;
else
    calc_c = 1;
end
    
%% main

if da == d
    L = (alpha - 1)' * log(X);
else
    L = (alpha - 1)' * sum(log(X), 1);
end

if calc_c
    c = mvbetaln(alpha, d);
end

if ~isequal(c, 0)
    if isscalar(c)
        L = L - c;
    else
        L = bsxfun(@minus, L, c.');
    end
end    


