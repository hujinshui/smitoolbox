function X = gsample(mu, sigma, n, i, rstream)
% Samples from (multivariate) Gaussian distributions
%
%   X = gsample(mu, sigma);
%       draws a sample from the Gaussian distribution whose mean and
%       covariance are respectively given by mu and sigma.
%
%       Here, mu is a d x m matrix (or a zero scalar to indicate zero
%       mean), and sigma is an object representing a covariance matrix
%       (such as the objects of classes udmat, dmat, gsymat, etc).
%       sigma.n can be 1 (the covariance is shared by all distributions),
%       or m (each distribution has its own covariance).
%
%       If m == 1, then X is a d x 1 column representing an obtained
%       sample. If m > 1, then X is a d x m matrix, comprised of m samples.
%       In particular, X(:,k) is from the k-th distribution.
%
%   X = gsample(mu, sigma, n);
%       If n is a scalar, then it draws n samples from each distribution 
%       given by mu and sigma.
%       
%       If n is a vector of length m, then it draws n(k) samples from
%       the k-th distribution.
%
%   X = gsample(mu, sigma, n, i);
%       If both n and i are scalars, then it draws n samples from the
%       i-th distribution.
%
%       The i can also be a vector. In this case, if n is a scalar, 
%       it draws n samples from each of the distribution selected by i.
%       n can also be a vector of the same size as i, then it draws
%       n(k) samples from the i(k)-th distribution.
%
%   X = gsample(mu, sigma, n, [], rstream);
%   X = gsample(mu, sigma, n, i, rstream);
%       one can further specify the random number stream for sampling.
%

%   History
%   -------
%       - Created by Dahua Lin, on Nov 16, 2010
%

%% verify input

if ~(isfloat(mu) && ndims(mu) == 2)
    error('gsample:invalidarg', 'mu should be a numeric matrix.');
end

if ~isobject(sigma)
    error('gsample:invalidarg', ...
        'sigma should be an object representing the covariance.');
end
sn = sigma.n;

if isequal(mu, 0)
    d = sigma.d;
    m = sn;
else
    [d, m] = size(mu);    
    if ~(sigma.d == d && (sn == 1 || sn == m))
        error('gsample:invalidarg', ...
            'The size of sigma is not consistent with that of mu.');
    end
end

if nargin < 3
    n = 1;
else
    if ~(isnumeric(n) && (isscalar(n) || isvector(n)))
        error('gsample:invalidarg', 'n should be either a scalar or a vector.');
    end
end
if nargin < 4; i = []; end    
if nargin < 5; rstream = []; end


%% main

% count number

if isempty(i)    
    if isscalar(n)
        N = m * n;
    else
        if numel(n) ~= m
            error('gsample:invalidarg', 'The size of n is invalid.');
        end
        N = sum(n);
    end

else % i not empty            
    m = numel(i);
    
    if isscalar(n)
        N = m * n;
    else
        if numel(n) ~= m
            error('gsample:invalidarg', 'The size of n is invalid.');
        end
        N = sum(n);
    end    
end

if m > 1
    if isempty(i)
        i = 1 : m;
    end
    if isscalar(n) && n ~= 1
        n = constmat(1, m, n);
    end
end 
    

% generate ~ N(0, 1) 

if isempty(rstream)
    X = randn(d, N);
else
    X = randn(rstream, d, N);
end

% linear transform with chol(C)

if sn == 1
    X = sigma.choltrans(X);
    
else  % sn = m > 1
    if isequal(n, 1)
        for k = 1 : m
            X(:,k) = sigma.take(i(k)).choltrans(X(:,k));
        end
    else
        b = 0;
        for k = 1 : m
            ci = b + (1:n(k));
            X(:, ci) = sigma.take(i(k)).choltrans(X(:, ci));
            b = b + n(k);
        end
    end
end

% add mu

if ~isequal(mu, 0)
    if isequal(n, 1)
        if isempty(i)
            X = X + mu;
        else
            X = X + mu(:, i);
        end
        
    elseif m == 1
        if isempty(i)
            X = bsxfun(@plus, X, mu);
        else
            X = bsxfun(@plus, X, mu(:,i));
        end
        
    else
        % i must not be empty when m > 1
        b = 0;
        for k = 1 : m
            ci = b + (1:n(k));
            X(:, ci) = bsxfun(@plus, X(:, ci), mu(:, i(k)));
            b = b + n(k);
        end
    end
end


