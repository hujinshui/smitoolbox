function X = gaussd_sample(G, n, Z)
% Samples from (multivariate) Gaussian distributions
%
%   X = gaussd_sample(G);
%   X = gaussd_sample(G, n);
%       Draws n samples from a Gaussian distribution G whose mean and 
%       covariance are respectively given by mu and C.
%
%       Input arguments:
%       - G:        a gaussd struct with G.n == 1.
%       - n:        the number of sameples to be acquired from the model.
%                   (If n is omitted, then n is assumed to be 1).
%
%       Outputs:
%       - X:        a d x n matrix comprised of the generated samples
%                   as columns.
%
%   X = gaussd_sample(G, [], Z);
%
%       Transforms the samples drawn from a standard Gaussian distribution
%       stored as columns of Z to the samples from G.
%
%       Z should be a d x n matrix.
%

%
%   History
%   -------
%       - Created by Dahua Lin, on Aug 17, 2011
%       - Modified by Dahua Lin, on Aug 25, 2011
%       - Modified by Dahua Lin, on Sep 27, 2011
%       - Modified by Dahua Lin, on Nov 30, 2011
%       - Modified by Dahua Lin, on Dec 5, 2011
%

%% verify input arguments

if ~(is_gaussd(G) && G.n == 1)
    error('gaussd_sample:invalidarg', ...
        'G should be a gaussd struct with G.n == 1.');
end

if nargin < 2
    n = 1;
end


%% main skeleton

d = G.d;

if ~isempty(n)
    Z = randn(d, n);
else
    if ~(isfloat(Z) && isreal(Z) && ndims(Z) == 2 && size(Z,1) == d)
        error('gaussd_sample:invalidarg', ...
            'Z should be a d x n matrix.');
    end
end
    
if G.ty == 'm'
    X = gsample_m(Z, G.mu, G.C);
else
    X = gsample_c(Z, G.h, G.J);
end

%% core functions

function X = gsample_m(X, mu, C)

[d, n] = size(X);
ty = C.ty;
v = C.v;

if ty == 's' || ty == 'd'   
    if ~isequal(v, 1);
        if isscalar(v) || n == 1
            X = X .* sqrt(v);
        else
            X = bsxfun(@times, X, sqrt(v));
        end
    end
        
elseif ty == 'f'
    
    L = chol(v, 'lower');
    X = L * X;
end

X = add_mu(d, n, X, mu);


function X = gsample_c(X, h, J)

[d, n] = size(X);
ty = J.ty;
v = J.v;

if ty == 's' || ty == 'd'
    
    if ~isequal(v, 1);
        v = 1 ./ v;
        if ~isequal(h, 0)
            mu = h .* v;
        end

        if isscalar(v) || n == 1
            X = X .* sqrt(v);
        else
            X = bsxfun(@times, X, sqrt(v));
        end
    end
        
elseif ty == 'f'
    L = chol(v);
    g = L' \ h;
    A = L \ [X g];
    X = A(:, 1:n);
    mu = A(:, n+1);
    
end

X = add_mu(d, n, X, mu);



%% auxiliary functions

function X = add_mu(d, n, X, mu)

if ~isequal(mu, 0)
    if d == 1 || n == 1
        X = mu + X;
    else
        X = bsxfun(@plus, X, mu);
    end
end

