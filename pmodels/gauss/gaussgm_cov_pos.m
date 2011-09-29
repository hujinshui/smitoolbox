function [S_x, nu_x, C] = gaussgm_cov_pos(cf, cpri, X, mu, Z, tied, op, skipverify)
% Infer the posterior of covariance matrix 
%
%   [S, nu] = gaussgm_cov_pos(cf, cpri, X, mu);
%   [S, nu] = gaussgm_cov_pos(cf, cpri, X, mu, Z);
%   [S, nu] = gaussgm_cov_pos(cf, cpri, X, mu, Z, tied);
%
%       Computes the parameters of the posterior distribution, according
%       to cf, the form of the covariance.
%
%       Specifically,
%
%       cf == 's':      cpri:   an object of invgammad with cpri.dim == 1
%                       S:      beta parameter of the posterior 
%                       nu:     alpha parameter of the posterior
%
%       cf == 'd':      cpri:   an object of invgammad with cpri.dim == d
%                       S:      beta parameters of the posterior
%                       nu:     alpha parameter of the posterior
%
%       cf == 'f':      cpri:   an object of invwishartd with cpri.dim == d
%                       S:      the inverse scalae matrix of the posterior
%                       nu:     the degree of freedom of the posterior
%
%       Other inputs:
%       - X:            the observation matrix [d x n]
%
%       - mu:           the component-wise mean vectors [d x K]
%
%       - Z:            the labeling, which can be 
%                       - empty, indicating no grouping or weighting
%                       - K x n weighting matrix
%                       - a cell array of grouped indices
%
%       - tied:         whether the covariance of different components
%                       are tied.
%
%   [S, nu, C] = gaussgm_cov_pos(cf, cpri, X, mu);
%   [S, nu, C] = gaussgm_cov_pos(cf, cpri, X, mu, Z);
%   [S, nu, C] = gaussgm_cov_pos(cf, cpri, X, mu, Z, tied);
%
%       These statements further compute the maximum a posteriori 
%       estimation of covariance, and return C in form of pdmat struct.
%
%   [S, nu, C] = gaussgm_cov_pos(cf, cpri, X, mu, Z, tied, 'sample');
%
%       Draws a sample from the posterior distribution and returns it
%       as C.
%

% Created by Dahua Lin, on Sep 28, 2011
%

%% verify inputs

if nargin < 5; Z = []; end
if nargin < 6; tied = false; end
if nargin < 7; op = []; end

if ~skipverify
    [d, n, K, zty, tied, samp] = verify_args(cf, cpri, X, mu, Z, tied, op);
else
    [d, n] = size(X);
    K = size(mu, 2);
    if isempty(Z)
        zty = 0;
    elseif isnumeric(Z)
        zty = 1;
    else
        zty = 2;
    end
    tied = tied || (K == 1);
    samp = (~isempty(op) && strcmpi(op, 'sample'));
end

if nargout < 3
    mapest = false;
    samp = false;
else
    mapest = ~samp;
end


%% main

% Compute S0 and nu0

if cf == 's' || cf == 'd'

    if zty == 0
        S_x = vscatter(X, mu, 1, 'v');
        nu_x = n;
    elseif zty == 1
        S_x = vscatter(X, mu, Z, 'v');
        nu_x = sum(Z, 2)';
    else
        S_x = zeros(d, K);
        nu_x = zeros(1, K);
        for k = 1 : K
            S_x(:,k) = vscatter(X, mu(:,k), Z{k}, 'v');
            nu_x(k) = numel(Z{k});
        end
    end    
    
    if cf == 's' && d > 1
        S_x = sum(S_x, 1) * (1 / d);
    end
    
    if tied && K > 1
        S_x = sum(S_x, 2);
        nu_x = sum(nu_x);
    end    
    
else
    
    if zty == 0
        S_x = vscatter(X, mu, 1, 'm');
        nu_x = n;
    elseif zty == 1
        S_x = vscatter(X, mu, Z, 'm');
        nu_x = sum(Z, 2)';
    else
        S_x = zeros(d, K);
        nu_x = zeros(1, K);
        for k = 1 : K
            S_x(:,k) = vscatter(X, mu(:,k), Z{k}, 'm');
            nu_x(k) = numel(Z{k});
        end
    end 
    
    if tied && K > 1
        S_x = sum(S_x, 3);
        nu_x = sum(nu_x);
    end    
    
end

% Incorporate prior with S and nu

if cf == 's' || cf == 'd'
    S_x = S_x * 0.5;
    nu_x = nu_x * 0.5;
    
    if isempty(cpri)
        S = S_x;
        nu = nu_x - 1;
        
    else
        if tied
            S = S_x + cpri.beta;
        else
            S = bsxfun(@plus, S_x, cpri.beta);
        end
        nu = nu_x + cpri.alpha;
    end
    
else    
    if isempty(cpri)
        S = S_x;
        nu = nu_x - (d + 1);        
    else
        S = S_x + pdmat_fullform(cpri.Phi);
        nu = nu_x + cpri.deg;        
    end    
end     

% MAP estimation

if mapest
    if cf == 's' || cf == 'd'        
        if tied || d == 1
            C = S .* (1 ./ (nu + 1));
        else
            C = bsxfun(@times, S, 1 ./ (nu + 1));
        end        
    else
        if tied || d == 1
            C = S .* (1 ./ (nu + d + 1));
        else
            C = bsxfun(@times, S, 1 ./ (nu + d + 1));
        end
    end 
    C = pdmat(cf, d, C);
end

% Sampling

if samp
    if cf == 's' || cf == 'd'
        J = gamma_sample(nu, S, 1);
        C = 1 ./ J;        
    else
        J = wishart_sample(S, nu, 1);
        C = inv(J);
        C = 0.5 * (C + C');
    end   
    C = pdmat(cf, d, C);
end


%% Argument verification

function [d, n, K, zty, tied, samp] = verify_args(cf, cpri, X, mu, Z, tied, op)

if ~(ischar(cf) && isscalar(cf))
    error('gaussgm_cov_pos:invalidarg', 'cf should be a char scalar.');
end

if ~(isfloat(X) && ndims(X) == 2 && isreal(X))
    error('gaussgm_cov_pos:invalidarg', 'X should be a real matrix.');
end
[d, n] = size(X);

if ~(isfloat(mu) && ndims(mu) == 2 && isreal(mu))
    error('gaussgm_cov_pos:invalidarg', 'mu should be a real matrix.');
end
if size(mu, 1) ~= d
    error('gaussgm_cov_pos:invalidarg', ...
        'The dimension of X and that of mu are inconsistent.');
end
K = size(mu, 2);

if cf == 's'
    if ~isempty(cpri)
        if ~(isa(cpri, 'scale_invchi2d') && cpri.num == 1 && cpri.dim == 1)
            error('gaussgm_cov_pos:invalidarg', 'cpri is invalid.');
        end
    end
    
elseif cf == 'd'
    if ~isempty(cpri)
        if ~(isa(cpri, 'scale_invchi2d') && cpri.num == 1 && cpri.dim == d)
            error('gaussgm_cov_pos:invalidarg', 'cpri is invalid.');
        end
    end
    
elseif cf == 'f'
    if ~isempty(cpri)
        if ~(isa(cpri, 'invwishartd') && cpri.num == 1 && cpri.dim == d)
            error('gaussgm_cov_pos:invalidarg', 'cpri is invalid.');
        end
    end
    
else
    error('gaussgm_cov_pos:invalidarg', 'The value of cf is invalid.');
end

if isempty(Z)
    if K > 1
        error('gaussgm_cov_pos:invalidarg', ...
            'Z must be explicitly specified when K > 1.');
    end
    zty = 0;
elseif isnumeric(Z)
    if ~(isfloat(Z) && isreal(Z) && isequal(size(Z), [K n]))
        error('gaussgm_cov_pos:invalidarg', ...
            'The weighting matrix must be a matrix of size K x n.');
    end
    zty = 1;
elseif iscell(Z)
    if numel(Z) ~= K
        error('gaussgm_cov_pos:invalidarg', ...
            'The number of cells in Z should be K.');
    end
    zty = 2;
else
    error('gaussgm_cov_pos:invalidarg', 'The form of Z is invalid.');
end

if ~(islogical(tied) && isscalar(tied))
    error('gaussgm_cov_pos:invalidarg', 'tied should be a logical scalar.');
end
tied = tied || (K == 1);

if isempty(op)
    samp = false;
else
    if ~(ischar(op) && strcmpi(op, 'sample'))
        error('gaussgm_cov_pos:invalidarg', ...
            'The last argument (op) is invalid.');
    end
    samp = true;
end
    
if samp && isempty(cpri)
    error('gaussgm_cov_pos:invalidarg', ...
        'The prior must be explicitly given for sampling.');
end




