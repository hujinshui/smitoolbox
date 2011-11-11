function X = gsample(mu, C, n, op)
% Samples from (multivariate) Gaussian distributions
%
%   X = gsample(mu, C, n);
%       Draws n samples from a Gaussian distribution whose mean and 
%       covariance are respectively given by mu and C.
%
%       Input arguments:
%       - mu:       the mean vector [d x 1]
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
%

%% verify input arguments

if ~(isfloat(mu) && isreal(mu) && ndims(mu) == 2 && size(mu, 2) == 1)
    error('gsample:invalidarg', ...
        'The 1st arg to gsample should be a floating-point real vector.');
end

if ~(is_pdmat(C))
    error('gsample:invalidarg', ...
        'The 2nd arg to gsample should be a pdmat struct.');
end

d = C.d;
if ~(size(mu, 1) == d || isequal(mu, 0))
    error('gsample:invalidarg', 'The dim of mu and C are inconsistent.');
end

if ~(isnumeric(n) && isscalar(n) && n == fix(n) && n >= 0)
    error('gsample:invalidarg', ...
        'n should be a numeric real integer scalar.');
end

use_ip = 0;
if nargin >= 4
    if ~(ischar(op) && strcmpi(op, 'ip'))
        error('gsample:invalidarg', ...
            'The 4th arg to gsample can only be ''ip'' if given.');
    end
    use_ip = 1;
end


%% main

if use_ip
    h = mu;
    J = C;
    
    C = pdmat_inv(J);    
    
    if isequal(h, 0)
        mu = 0;
    else
        mu = pdmat_mvmul(C, h);
    end    
end


X = pdmat_choltrans(C, randn(d, n));

if ~isequal(mu, 0)
    X = bsxfun(@plus, X, mu);
end






