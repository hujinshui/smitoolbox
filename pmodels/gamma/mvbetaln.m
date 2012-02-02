function v = mvbetaln(alpha, d)
%MVBETALN Logarithm of multivariate beta function
%
%   A multivariate beta function is defined to be
%
%       f(a) = prod(gamma(a)) / gamma(sum(a)).
%
%
%   v = mvbetaln(alpha);
%
%       Evaluates the logarithm of multivariate beta function.
%       
%       Here, alpha should be a matrix of size d x n, and in the output,
%       the size of v is 1 x n. Particularly, v(i) corresponds to 
%       the parameter given by alpha(:,i).
%
%   v = mvbetaln(alpha, d);
%
%       Evaluates the logarithm of (symmetric) multivariate beta function.
%
%       When alpha is a row vector of size 1 x n, this is equivalent to
%       mvbetaln(repmat(alpha, d, 1)).
%
%       When alpha is a matrix of size d x n, this is equivalent to 
%       mvbetaln(alpha, d).
%

% Created by Dahua Lin, on Dec 26, 2011
%

%% verify input

if ~(isfloat(alpha) && isreal(alpha) && ndims(alpha) == 2)
    error('mvbetaln:invalidarg', 'alpha should be a real matrix.');
end
da = size(alpha, 1);

if nargin < 2
    d = da;
else
    if ~(isscalar(d) && isnumeric(d))
        error('mvbetaln:invalidarg', 'd should be a numeric scalar.');
    end    
    if ~(da == 1 || da == d)
        error('mvbetaln:invalidarg', 'd is inconsistent with alpha.');
    end
end

%% main

if d == 1
    v = zeros(1, size(alpha, 2));
else
    if da == 1
        v = gammaln(alpha) * d - gammaln(alpha * d);
    else
        v = sum(gammaln(alpha), 1) - gammaln(sum(alpha, 1));
    end
end

