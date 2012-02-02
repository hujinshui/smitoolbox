function G = gaussd_mle(X, W, cform, tie_cov)
% Performs Maximum likelihood estimation of Gaussian distributions
%
%   G = gaussd_mle(X);
%   G = gaussd_mle(X, W);
%   G = gaussd_mle(X, W, cform);
%   G = gaussd_mle(X, W, cform, tie_cov);
%
%       performs maximum likelihood estimation of Gaussian models based on
%       a (weighted) set of samples given by columns of X.
%
%       Input arguments:
%       - X:        the sample matrix of size d x n. Each column is a
%                   sample. 
%
%       - W:        The sample weights. 
%                   It can be omitted, empty, or a K x n matrix.
%                   If omitted or empty, then all samples are assumed to 
%                   have the same weight. If W is a K x n matrix, then
%                   K distributions are to be estimated, and the k-th one
%                   is estimated based on the weights given in W(k,:).
%
%       - cform:    the char indicating the form of covariance matrix.
%                   It can take either of the following values:
%                   - 's':  isotropic covariance in form of c * I
%                   - 'd':  diagonal covariance
%                   - 'f':  full covariance form
%                   If omitted, cform is assumed to be 'f'.
%
%       - tie_cov:  If multiple distributions are to be estimated, whether
%                   their covariance is tied to the same one.
%                   If omitted, tie_cov is assumed to be false.
%       

%   History
%   -------
%       - Created by Dahua Lin, on Nov 14, 2010
%       - Modified by Dahua Lin, on Aug 25, 2010
%       - Modified by Dahua Lin, on Dec 6, 2010
%

%% verify input
            
if isfloat(X) && ndims(X) == 2
    [d, n] = size(X);
else
    error('gaussd_mle:invalidarg', 'X should be a d x n numeric matrix.');
end

if nargin < 2 || isempty(W)
    W = [];
    K = 1;
else
    if ~(isfloat(W) && isreal(W) && ...
            (isscalar(W) || (ndims(W)==2 && size(W,2) == n)))
        error('gaussd_mle:invalidarg', ...
            'W should be a scalar or a real matrix with n columns.');
    end
    K = size(W, 1);
end

if nargin < 3
    cform = 'f';
else
    if ~(ischar(cform) && isscalar(cform) && any(cform == 'sdf'))
        error('gaussd_mle:invalidarg', ...
            'cform should be either ''s'', ''d'', or ''f''.');
    end
end

if nargin < 4
    tie_cov = false;
else
    if ~(islogical(tie_cov) && isscalar(tie_cov))
        error('gaussd_mle:invalidarg', 'tie_cov should be a logical scalar.');
    end
end

%% main

% preparation

if isempty(W)
    Wt = [];
else
    % normalize the weights
    Wt = W.';    
    
    sw = sum(Wt, 1);
    if issparse(sw)
        sw = full(sw);
    end
    Wt = bsxfun(@times, Wt, 1 ./ sw);
    sw = sw.' / sum(sw);
end

% estimate mean vectors

mu = mean_w(X, Wt);

% estimate variance / covariance

switch cform
    case 's'
        if d == 1
            ex2 = mean_w(X .* X, Wt);
        else
            ex2 = mean_w(dot(X, X, 1), Wt);
        end
        v = ex2 - dot(mu, mu, 1);
        if K > 1 && tie_cov
            v = mean_w(v, sw);
        end
        if d == 1
            C = pdmat('s', 1, v);
        else
            C = pdmat('s', d, v * (1/d));
        end
        
    case 'd'
        ex2 = mean_w(X.^2, Wt);
        v = ex2 - mu .^ 2;
        if K > 1 && tie_cov
            v = mean_w(v, sw);
        end
        C = pdmat('d', d, v);
                
    case 'f'
        if K == 1
            C = calc_cov(X, mu, Wt);
        else
            if tie_cov
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
        C = pdmat('f', d, C);
end

% generate Gaussian model

G = gaussd('m', mu, C);


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
    if ~issparse(Wt)
        Exx = X * bsxfun(@times, X', Wt);
    else
        [I, ~, w] = find(Wt);
        X = X(:, I);
        Exx = X * bsxfun(@times, X', w);
    end
end

C = Exx - mu * mu';
C = (C + C') * 0.5;

