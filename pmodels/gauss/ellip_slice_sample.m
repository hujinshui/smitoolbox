function X = ellip_slice_sample(x0, llikfun, V, intv, n)
%ELLIP_SLICE_SAMPLE Elliptical Slice Sampling
%
%   x = ELLIP_SLICE_SAMPLE(x0, llikfun, v);
%  
%       Draws the next sample from the posterior formulated as below
%
%           x ~ N(x | 0, Sigma) * llik( x );
%
%       Here, Sigma is the prior covariance matrix, and llik is an
%       aribtrary likelihood function, which is given.
%
%       Input arguments:
%       - x0:           The current sample
%
%       - v:            A sample vector drawn from N(0, Sigma).
%
%       - llikfun:      The log-likelihood function, which evaluates the 
%                       log-likelihood of a given parameter. 
%
%       This statement outputs a sample obtained by one iteration of 
%       the sampling procedure.
%
%   X = ELLIP_SLICE_SAMPLE(x0, V, llikfun, intv, n);
%
%       Generates n samples by running the elliptical slice sampling
%       algorithm. The function runs intv iterations before generating
%       each sample.
%
%       Note that the size of V, the matrix comprised of all pre-sampled
%       vectors from the prior, should be d x (intv x n), meaning there
%       are intv x n columns in V.
%

% Created by Dahua Lin, on Feb 25, 2012
%

%% verify input arguments

if ~(isfloat(x0) && isreal(x0) && ndims(x0) == 2 && size(x0,2) == 1)
    error('ellip_slice_sample:invalidarg', ...
        'x0 should be a real vector.');
end
d = size(x0, 1);

if ~(isa(llikfun, 'function_handle'))
    error('ellip_slice_sample:invalidarg', ...
        'llikfun should be a function handle.');
end

if ~(isfloat(V) && isreal(V) && ismatrix(V) && size(V,1) == d)
    error('ellip_slice_sample:invalidarg', ...
        'V should be a real vector or matrix with size(V,1) == d.');
end

if nargin < 4
    intv = 1;
else
    if ~(isnumeric(intv) && isscalar(intv) && intv >= 1)
        error('ellipse_slice_sample:invalidarg', ...
            'intv should be a positive number.');
    end
end

if nargin < 5
    n = 1;
else
    if ~(isnumeric(n) && isscalar(n) && n >= 1)
        error('ellipse_slice_sample:invalidarg', ...
            'n should be a positive number.');
    end
end

if size(V, 2) ~= intv * n
    error('ellip_slice_sample:invalidarg', ...
        'The number of columns in V is incorrect.');
end
    

%% main

x = x0;
likv = llikfun(x);

if n == 1
    if intv == 1
        X = run_ess(x, likv, V, llikfun);
    else
        for j = 1 : intv
            [x, likv] = run_ess(x, likv, V(:,j), llikfun);
        end
        X = x;
    end
    
else
    X = zeros(d, n);    
    if intv == 1        
        for i = 1 : n
            [x, likv] = run_ess(x, likv, V(:,i), llikfun);
            X(:,i) = x;
        end
    else        
        for i = 1 : n
            for j = 1 : intv
                v = V(:, (i-1)*intv+j);
                [x, likv] = run_ess(x, likv, v, llikfun);
            end
            X(:,i) = x;
        end
    end
end


%% core function

function [x, likv] = run_ess(x, likv, v, llikfun)

% log-lik thres

logy = likv + log(rand());

% initialize

t = (2*pi) * rand();
rb = t;
lb = rb - (2*pi);

tx = x * cos(t) + v * sin(t);
tlik = llikfun(tx);

% loop

while tlik < logy
    
    % shrink the range
    if t < 0
        lb = t;
    else
        rb = t;
    end
    
    % pick a new t
    
    t = lb + rand() * (rb - lb);
    tx = x * cos(t) + v * sin(t);  
    tlik = llikfun(tx);
end

x = tx;
likv = tlik;

