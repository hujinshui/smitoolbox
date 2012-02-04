function M = ppca_mle(X, w, q, varargin)
% Performs Maximum Likelihood estimation of PPCA from data
%
%   M = ppca_mle(X, [], q, ...);
%   M = ppca_mle(X, w, q, ...);
%
%       performs Maximum Likelihood estimation of the PPCA
%       model from data.
%
%       Input arguments:
%       - X:    The data matrix of size d x n, of which each
%               column is s sample.
%       - w:    The sample weights, a vector of length n.
%               If all samples have the same weight, it can be empty.
%       - q:    the dimension of the latent space. It should
%               have q < min(d, n).
%
%       One can specify further options to control the
%       estimation, in form of name/value pairs.
%
%       - 'method':     The method used to do the training,
%                       which can be
%                       - 'cov':    by computing the covariance
%                                   matrix, and compute the
%                                   eigenvectors of it.
%                       - 'svd':    by doing SVD, this can be
%                                   faster when n < d.
%                       The default is 'cov'.
%
%       - 'mu':         Pre-computed mean vector, which can be 
%                       either of the following:
%                       - []:   no pre-computed vector. The function
%                               will compute mu.
%                       - a d x 1 mean vector.
%                       - a zero scalar, indicating zero mean.
%

%   History
%   -------
%       - Created by Dahua Lin, on Nov 20, 2010
%       - Modified by Dahua Lin, on Nov 3, 2011
%       - Modified by Dahua Lin, on Dec 27, 2011
%


%% verify arguments

if ~(isfloat(X) && ndims(X) == 2 && isreal(X))
    error('probpca:mle:invalidarg', ...
        'X should be a real matrix.');
end
[d, n] = size(X);

if ~isempty(w)
    if ~(isfloat(w) && isreal(w) && isvector(w) && numel(w) == n)
        error('probpca:mle:invalidarg', ...
            'w should be a vector of length n.');
    end
    if size(w, 2) > 1
        w = w.';
    end
end

if ~(isnumeric(q) && isscalar(q) && q == fix(q) && ...
        q >= 1 && q < d && q < n);
    error('probpca:mle:invalidarg', ...
        'q should be a positive integer less than min(d, n).');
end

[method, mu] = chk_opts(d, varargin);

%% main

% mean vector

if isempty(mu)
    if isempty(w)
        mu = sum(X, 2) * (1 / n);
    else
        sw = sum(w);
        mu = (X * w) * (1 / sw);
    end    
end

if isequal(mu, 0)
    Z = X;
else
    Z = bsxfun(@minus, X, mu);
end

% basis and eigenvalues

switch method
    case 'cov'
        if isempty(w)
            C = (Z * Z') * (1/n);
        else
            C = (Z * bsxfun(@times, Z', w)) * (1/sw);
            C = 0.5 * (C + C');
        end
        [U, evs] = eig(C);
        evs = diag(evs);
        
    case 'svd'
        if isempty(w)
            [U, svs] = svd(Z, 0);
            svs = diag(svs);
            evs = (svs.^2) * (1/n);
        else
            [U, svs] = svd(bsxfun(@times, Z, sqrt(w)'), 0);
            svs = diag(svs);
            evs = (svs.^2) * (1/sw);
        end
        
    otherwise
        error('ppca_mle:invalidarg', 'The method %s is invalid.', method);
end

% make struct

d = size(U, 1);
[evs, si] = sort(evs, 1, 'descend');
qevs = evs(1:q);
qevs = max(qevs, 1e-12);
B = U(:, si(1:q));

se2 = sum(evs(q+1:end)) / (d-q);
se2 = max(se2, 1e-12);

M = ppca_model(B, sqrt(qevs - se2), sqrt(se2), mu);
            
            
%% sub functions

function [method, mu] = chk_opts(d, params)

method = 'cov';
mu = [];

if ~isempty(params)
    onames = params(1:2:end);
    ovals = params(2:2:end);
    
    if ~(numel(onames) == numel(ovals) && iscellstr(onames))
        error('ppca_mle:invalidarg', ...
            'The name/value list for options is invalid.');
    end
    
    for i = 1 : numel(onames)
        cn = onames{i};
        cv = ovals{i};
        switch lower(cn)
            case 'method'
                if ~(ischar(cv) && ...
                        (strcmp(cv, 'cov') || strcmp(cv, 'svd')))
                    error('ppca_mle:invalidarg', ...
                        'The method should be either ''cov'' or ''svd''.');
                end
                method = cv;
            case 'mu'
                if ~isempty(cv)
                    if ~( isequal(cv, 0) || (isfloat(cv) && isreal(cv) && ...
                            isequal(size(cv), [d 1])) )
                        error('ppca_mle:invalidarg', ...
                            'mu should be either 0 or a d x 1 real vector.');
                    end
                    mu = cv;
                end                    
                    
            otherwise
                error('ppca_mle:invalidarg', ...
                    'Unknown option name %s', cn);
        end
    end
end

