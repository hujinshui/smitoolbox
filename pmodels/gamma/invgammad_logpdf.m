function L = invgammad_logpdf(alpha, beta, X, c)
%invgammad_logpdf Evaluate log probability density of inverse gamma distribution
%
%   log f(x; alpha, beta) = 
%       - (alpha + 1) * log(x) - beta / x + const
%
%   with const = alpha * log(beta) - gammaln(alpha)
%
%   L = invgammad_logpdf(alpha, beta, X);
%   L = invgammad_logpdf(alpha, beta, X, 0);
%   L = invgammad_logpdf(alpha, beta, X, c);
%
%       evaluates the log probability density function of inverse gamma
%       distribution(s) at given samples.
%
%       Inputs:
%       - alpha:    the shape parameter(s)
%       - beta:     the scale parameter(s)
%       - X:        the sample matrix [d x n]
%
%       - c:        the constant term in the log-pdf. 
%
%                   If c is not provided if will be evaluated in the 
%                   function. When this function is invoked multiple times 
%                   with the same set of distributions, it is advisable 
%                   to pre-compute it using gammad_const.
%
%                   One can also set c to 0, which means only computing
%                   the linear part without adding the constant.
%
%       When d > 1, this indicates a multi-dimensional inverse gamma 
%       distribution with independent components.
%
%       alpha and beta can respectively be either of the following:
%       a scalar, d x 1 vector, 1 x m vector, or d x m matrix. 
%       When m > 1, it indicates there are multiple distributions, 
%       the log-pdf values with respect to all distributions are evaluated
%       for each sample.
%
%       The sizes of alpha and beta need not be the same, but they have 
%       to be compatible with each other in bsxfun sense.
%
%       Outputs:
%       - L:        the log-pdf value matrix, of size m x n.
%                   Particularly, L(k, i) is the log-pdf at the i-th sample
%                   with respect to the k-th distribution.
%
%
%

% Created by Dahua Lin, on Dec 26, 2011
%

%% verify input arguments

if ~(isfloat(X) && isreal(X) && ndims(X) == 2)
    error('invgammad_logpdf:invalidarg', 'X should be a real matrix.');
end
if ~(isfloat(alpha) && isreal(alpha) && ndims(alpha) == 2)
    error('invgammad_logpdf:invalidarg', 'alpha should be a real matrix.');
end
if ~(isfloat(beta) && isreal(beta) && ndims(beta) == 2)
    error('invgammad_logpdf:invalidarg', 'beta should be a real matrix.');
end

dx = size(X, 1);
[da, ma] = size(alpha);
[db, mb] = size(beta);

if ~( (da == 1 || da == dx) && (db == 1 || db == dx) && ...
        (ma == 1 || mb == 1 || ma == mb))
    error('invgammad_logpdf:invalidarg', 'The size of alpha or beta is invalid.');
end

m = max(ma, mb);

if nargin >= 4
    if ~( isequal(c, 0) || ...
            (isfloat(c) && isreal(c) && isequal(size(c), [1 m])) )
        error('invgammad_logpdf:invalidarg', ...
            'c should be either zero or a 1 x m real vector.');
    end
    calc_c = 0;
else
    calc_c = 1;
end


%% Evaluate

% first term: (alpha - 1) log(x)

if da == dx
    T1 = (alpha + 1)' * log(X);
else
    T1 = (alpha + 1)' * sum(log(X), 1);
end

% second term: x / beta

if db == dx
    T2 = beta' * (1 ./ X);
else
    T2 = beta' * sum(1 ./ X, 1);
end

% combine terms

if size(T1, 1) == size(T2, 1)
    L = - (T1 + T2);
else
    L = - bsxfun(@plus, T1, T2);
end

% add constants

if calc_c
    c = invgammad_const(alpha, beta);
    if da < dx && db < dx
        c = c * dx;
    end
end

if ~isequal(c, 0)
    if m == 1
        L = L + c;
    else
        L = bsxfun(@plus, L, c.');
    end    
end
    

