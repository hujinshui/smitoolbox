function X = gsample(mu, C, n, op)
% Samples from (multivariate) Gaussian distributions
%
%   X = gsample(mu, C, n);
%       Draws n samples from a Gaussian distribution whose mean and 
%       covariance are respectively given by mu and C.
%
%       Input arguments:
%       - mu:       the mean vector [d x 1], or just a zero.
%       - C:        the covariance given in pdmat struct.
%       - n:        the number of sameples to be acquired from the model.
%   
%   X = gsample(h, J, n, 'ip');
%       Performs sampling based on information parameters
%

%
%   History
%   -------
%       - Created by Dahua Lin, on Aug 17, 2011
%       - Modified by Dahua Lin, on Aug 25, 2011
%       - Modified by Dahua Lin, on Sep 27, 2011
%       - Modified by Dahua Lin, on Nov 30, 2011
%

%% verify input arguments

if isfloat(C)
    if ~(isreal(C) && size(C,1) == size(C,2))
        error('gsample:invalidarg', 'C or J should be a real square matrix.');
    end
    if isscalar(C)
        ty = 's';
        d = size(mu, 1);
    else
        ty = 'f';
        d = size(C, 1);
    end
    v = C;
    
elseif is_pdmat(C)
    if C.n ~= 1
        error('gsample:invalidarg', 'C.n or J.n must equal 1 for gsample.');
    end
    ty = C.ty;
    v = C.v;
    d = C.d;

else
    error('gsample:invalidarg', 'The 2nd arg to gsample is invalid.');
end

use_ip = 0;
if nargin >= 4
    if ~strcmp(op, 'ip')
        error('gsample:invalidarg', ...
            'The 4th arg to gsample can only be ''ip'' if given.');
    end
    use_ip = 1;
end


%% main

X = randn(d, n);

if ty == 's' || ty == 'd'   
    if ~isequal(v, 1);
        if use_ip
            v = 1 ./ v;
            if ~isequal(mu, 0)
                mu = mu .* v;
            end
        end
        if isscalar(v) || n == 1
            X = X .* sqrt(v);
        else
            X = bsxfun(@times, X, sqrt(v));
        end
    end
        
elseif ty == 'f'
    
    if ~use_ip
        L = chol(v, 'lower');
        X = L * X;        
    else
        L = chol(v);
        g = L' \ mu;
        A = L \ [X g];
        X = A(:, 1:n);
        mu = A(:, n+1);
    end    
    
end

X = add_mu(d, n, X, mu);




%% sub functions

function X = add_mu(d, n, X, mu)

if ~isequal(mu, 0)
    if d == 1 || n == 1
        X = mu + X;
    else
        X = bsxfun(@plus, X, mu);
    end
end

