function X = gsample(cf, mu, C, n)
% Samples from (multivariate) Gaussian distributions
%
%   X = gsample(cf, mu, C, n);
%       Draws n samples from a Gaussian distribution whose mean and 
%       covariance are respectively given by mu and C.
%
%       Input arguments:
%       - cf:       the form of covariance
%                   which can be either 's', 'd', or 'f'.
%       - mu:       the mean vector [d x 1]
%       - C:        the covariance given in the specified form:
%                   - cf == 's':    C is a scalar that represents a
%                                   covariance matrix as C * eye(d).
%                   - cf == 'd':    C is a d x 1 vector that represents a
%                                   covariance matrix as diag(C).
%                   - cf == 'f':    C is a full covariance matrix.
%
%       - n:        the number of sameples to be acquired from the model.

% Created by Dahua Lin, on Aug 17, 2011
%

%% verify input arguments

if ~(ischar(cf) && isscalar(cf))
    error('gsample:invalidarg', ...
        'cf should be a char scalar.');
end

if ~(isfloat(mu) && isreal(mu) && ndims(mu) == 2 && size(mu, 2) == 1)
    error('gsample:invalidarg', ...
        'mu should be a floating-point real vector.');
end
d = size(mu, 1);

if ~(isfloat(C) && isreal(C) && ndims(C) == 2)
    error('gsample:invalidarg', ...
        'C should be a floating-point real matrix.');
end

if ~(isnumeric(n) && isscalar(n) && n == fix(n) && n >= 0)
    error('gsample:invalidarg', ...
        'n should be a numeric real integer scalar.');
end

%% main

if cf == 's' || cf == 'd' || cf == 'f'
    X = gmat_choltrans(cf, C, randn(d, n));
        
else
    error('gsample:invalidarg', ...
        'Invalid char for covariance form.');
end

if ~all(mu == 0)
    X = bsxfun(@plus, X, mu);
end

