function X = dird_sample(alpha, n)
%DIRD_SAMPLE Samples from a Dirichlet distribution
%
%   X = dird_sample(alpha);
%   X = dird_sample(alpha, n);
%   X = dird_sample(alpha, [d, n]);
%
%       Draws n samples from a Dirichlet distribution over the 
%       (d-1)-dimensional probability simplex, whose parameter
%       is given by alpha.
%
%       Input arguments:
%       - alpha:    the parameter of Dirichlet distribution. It is
%                   either a d x 1 column vector, or a scalar (for
%                   symmetric Dirichlet)
%
%       - d:        the sample dimension. Each sample is a vector
%                   of length d that sums to 1.
%
%       - n:        the number of samples to draw.
%                   

% Created by Dahua Lin, on Sep 27, 2011
%

%% verify input arguments

if ~( isfloat(alpha) && isreal(alpha) && ...
        (isscalar(alpha) || (ndims(alpha) == 2 && size(alpha,2) == 1)) )
    error('dird_sample:invalidarg', ...
        'alpha should be a real scalar or a real column vector.');
end
da = size(alpha, 1);

if nargin < 2
    n = 1;
    d = da;   
else    
    if ~(isnumeric(n) && (isscalar(n) || numel(n) == 2))
        error('dird_sample:invalidarg', ...
            'n should be a non-negative integer scalar or a numeric pair.');
    end
        
    if isscalar(n)
        d = da;
    else
        d = n(1);
        n = n(2);
        
        if ~(da == 1 || da == d)
            error('dird_sample:invalidarg', 'Inconsistent dimensions.');
        end
    end
end

%% main

if d > 1

    if isscalar(alpha)
        X = randg(alpha, d, n);
    else
        if n == 1
            X = randg(alpha);
        else
            X = randg(alpha(:, ones(1,n)));
        end
    end
    X = bsxfun(@times, X, 1 ./ sum(X, 1));
    
else % d == 1
    X = ones(1, n); 
    
end

