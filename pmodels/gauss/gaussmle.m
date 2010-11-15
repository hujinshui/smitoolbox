function G = gaussmle(X, W, varargin)
% Performs Maximum likelihood estimation of Gaussian models
%
%   G = gaussmle(X, W, ...);
%       performs maximum likelihood estimation of Gaussian models based on
%       the samples given by columns of X.
%
%       The samples can be weighted. The weights are given by rows of W.
%       Suppose there are n samples, then W can be a matrix of size 
%       K x n. Then the output G will be a gaussd object comprised of
%       K components. 
%
%       If the samples are not weighted, then one can set W to a scalar 1.
%
%       One can specify additional parameters
%
%       - 'cov_type':   the covariance type, which can be either of 
%                       the following: 
%                       - 'full':   full covariance matrix
%                       - 'diag':   diagonal covariance
%                       - 'iso':    isotropic covariance, i.e. eye(d) * v
%                       default = 'full'
%
%       - 'tie_cov':    whether to tie the covariance of all models
%                       together in estimation.
%                       default = false;
%
%       - 'use_ip':     whether to use information parameters
%                       default = false;
%

% Created by Dahua Lin, on Nov 14, 2010
%

%% verify input
            
if isfloat(X) && ndims(X) == 2
    [d, n] = size(X);
else
    error('gaussgm:estimate_mle:invalidarg', ...
        'X should be a d x n numeric matrix.');
end

if nargin < 2
    W = 1;
else
    if ~(isfloat(W) && ...
            (isscalar(W) || (ndims(W)==2 && size(W,2) == n)))
        error('gaussgm:estimate_mle:invalidarg', ...
            'W should be a scalar or a numeric matrix with n columns.');
    end
end

opts = struct('cov_type', 'full', 'tie_cov', false, 'use_ip', false);
if nargin > 2
    opts = parlist(opts, varargin{:});
end

%% preparation

switch opts.cov_type
    case 'iso'
        ctype = 1;
    case 'diag'
        ctype = 2;
    case 'full'
        ctype = 3;
    otherwise
        error('Invalid covariance type %s', opts.cov_type);
end        

if isscalar(W)
    Wt = [];
    K = 1;
else
    % normalize the weights
    sw = sum(W, 2);
    Wt = bsxfun(@times, W, 1 ./ sw).';
    sw = sw / sum(sw);
    K = size(W, 1);
end

%% do computation

% estimate mean vectors

mu = mean_w(X, Wt);

% estimate variance / covariance

if ctype == 1
    ex2 = mean_w(dot(X, X, 1), Wt);
    v = ex2 - dot(mu, mu, 1);    
    if K > 1 && opts.tie_cov
        v = mean_w(v, sw);
    end    
    sigma = udmat(d, v / d);
    
elseif ctype == 2
    ex2 = mean_w(X.^2, Wt);
    v = ex2 - mu .^ 2;
    if K > 1 && opts.tie_cov
        v = mean_w(v, sw);
    end
    sigma = dmat(v);
    
elseif ctype == 3
    if K == 1
        C = calc_cov(X, mu, Wt);
    else        
        if opts.tie_cov
            C = zeros(d, d);
            for k = 1 : K
                C = C + sw(k) * calc_cov(X, mu(:,k), Wt(:,k));
            end
        else
            C = zeros(d, d, K);
            for k = 1 : K
                C(:,:,k) = calc_cov(X, mu(:,k), Wt(:,k));
            end
        end
    end
    sigma = gsymat(C);
    
end

if opts.use_ip
    G = gaussd.from_mp(mu, sigma, 'ip');
else
    G = gaussd.from_mp(mu, sigma);
end


%% Auxiliary functions

function y = mean_w(x, w)

if isempty(w)
    y = sum(x, 2) * (1 / size(x,2));
else
    y = x * w;
end


function C = calc_cov(X, mu, Wt)

if isempty(Wt)
    Exx = (X * X') * (1 / size(X,2));
else
    Exx = X * bsxfun(@times, X', Wt);
end

C = Exx - mu * mu';
C = (C + C') * 0.5;

