function C = dird_cov(alpha, d)
%DIRD_COV the covariance matrix of a Dirichlet distribution
%
%   C = DIRD_COV(alpha);
%   C = DIRD_COV(alpha, d);
%
%       Computes the covariance matrix of a dirichlet distribution.
%       Here, alpha is either a scalar, or a vector of size d x 1.
%

% Created by Dahua Lin, on Dec 26, 2011
%

%% verify input arguments

if ~( isfloat(alpha) && isreal(alpha) && ...
        (isscalar(alpha) || (ndims(alpha) == 2 && size(alpha, 2) == 1)) )
    error('dird_cov:invalidarg', 'alpha should be a scalar or a column vector.');
end
da = size(alpha, 1);

if nargin < 2
    d = da;
else
    if ~(isscalar(d) && isnumeric(d))
        error('dird_cov:invalidarg', 'd should be a numeric scalar.');
    end    
    if ~(da == 1 || da == d)
        error('dird_cov:invalidarg', 'd is inconsistent with alpha.');
    end
end

%% main

if d == 1
    C = 0;
    return;
end

if da == 1
    a0 = alpha * d;
    
    s = 1 ./ (a0^2 * (a0 + 1));
    c1 = alpha .* (a0 - alpha) * s;
    c2 = - alpha^2 * s;
    
    C = constmat(d, d, c2);
    C(1:(d+1):d^2) = c1;
else
    a0 = sum(alpha, 1);
    
    s = 1 ./ (a0^2 * (a0 + 1));
    C = (alpha * alpha') * (-s);
    
    c1 = alpha .* (a0 - alpha) * s;
    C(1:(d+1):d^2) = c1;
end
    

    
