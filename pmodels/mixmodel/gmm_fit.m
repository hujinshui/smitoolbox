function [R, state] = gmm_fit(X, w, K, varargin)
%GMM_FIT Fits a Gaussian Mixture Model
%
%   R = gmm_fit(X, [], K, ...);
%   R = gmm_fit(X, w, K, ...);
%
%   [R, state] = gmm_fit( ... );
%
%       Fits a Gaussian mixture model to a (weighted) set of data.
%
%       Suppose there are n samples on a d dimensional space.
%
%       Input arguments:
%       - X:        the data matrix [d x n]
%       - w:        the sample weights, which can be either empty or
%                   a row vector of 1 x n.
%       - K:        the number of mixture components (K >= 2).      
%
%       Output arguments:
%       - R:        a struct with the following fields:
%                   - K:    the number of mixture components
%                   - Pi:   the prior distribution over mixture components
%                   - G:    the gaussd struct with G.d == d and G.n == K
%                   - Q:    the soft assignment matrix [K x n]
%
%       - state:    the fmm_std state object
%
%       One can further specify the following options in form of name
%       value list:
%
%       - 'cov_form':   the form of covariance:
%                           's':    isotropic covariance
%                           'd':    diagonal covariance
%                           'f':    full form covariance
%                       (default = 'f')
%
%       - 'initL':      the initial assignment, given by a label vector
%                       of size 1 x n. (default = [], indicating that
%                       a random initialization will be used).
%
%       - 'maxiters':   the maximum number of iterations in E-M 
%                       (default = 100)
%
%       - 'tol':        the tolerance of objective change upon convergence
%                       (default = 1e-6)
%
%       - 'display':    the level of displaying:
%                           'off':      no display
%                           'final':    display at the end
%                           'stage':    display per stage
%                           'eval':     display per objective evaluation
%                           'iter':     display per iteration
%                       (default = 'iter')
%
%       - 'tied_cov':   whether to tie covariance across components
%                       (default = false)
%
%       - 'pricount':   the prior count of each component (default = 0)
%
%

% Created by Dahua Lin, on Dec 27, 2011
%

%% verify input arguments

if ~(isfloat(X) && isreal(X) && ndims(X) == 2)
    error('gmm_fit:invalidarg', 'X should be a real matrix.');
end
[d, n] = size(X);

if ~isempty(w)
    if ~(isfloat(w) && isreal(w) && isvector(w) && numel(w) == n)
        error('gmm_fit:invalidarg', ...
            'w should be either empty or a real vector of length n.');
    end
    if size(w, 2) > 1
        w = w.';
    end
else
    w = [];
end

if ~(isnumeric(K) && isscalar(K) && K == fix(K) && K >= 2)
    error('gmm_fit:invalidarg', ...
        'K should be a positive integer with K >= 2.');
end

opts = parse_options(n, K, varargin);


%% main

if ~opts.tied_cov
    gm = gaussgm(d, opts.cov_form);
else
    gm = gaussgm(d, opts.cov_form, 'tied-cov');
end

gs = fmm_std('em', gm, [], opts.pricount);

if isempty(opts.initL)
    gs = gs.initialize(X, w, 'rand', K);
else
    gs = gs.initialize(X, w, 'labels', opts.initL);
end

vdopts = varinfer_options([], ...
    'maxiters', opts.maxiters, 'tol', opts.tol, 'display', opts.display);
ret = varinfer_drive(gs, vdopts);

% output

sol = ret.sol;

R.K = sol.K;
R.Pi = sol.Pi;
R.G = sol.params;
R.Q = sol.Q;

state = ret.state;



%% option parsing

function opts = parse_options(n, K, nvlist)

opts.cov_form = 'f';
opts.initL = [];
opts.maxiters = 100;
opts.tol = 1e-6;
opts.pricount = 0;
opts.tied_cov = false;
opts.display = 'off';


if ~isempty(nvlist)
    ns = nvlist(1:2:end);
    vs = nvlist(2:2:end);
    
    if ~(numel(ns) == numel(vs) && iscellstr(ns))
        error('gmm_fit:invalidarg', 'Invalid name value list.');
    end
    
    for i = 1 : numel(ns)
        
        cn = ns{i};
        cv = vs{i};
        
        switch lower(cn) 
            case 'cov_form'
                if ~(ischar(cv) && isscalar(cv) && any(cv == 'sdf'))
                    error('gmm_fit:invalidarg', ...
                        'cov_form should be either ''s'', ''d'', or ''f''.');
                end
                opts.cov_form = cv;
            
            case 'initl'
                if ~isempty(cv)
                    if ~(isnumeric(cv) && isreal(cv) && isequal(size(cv), [1, n]))
                        error('gmm_fit:invalidarg', ...
                            'initL should be a numeric vector of size 1 x n.');
                    end
                end
                opts.initL = cv;
            
            case 'maxiters'
                if ~(isnumeric(cv) && isscalar(cv) && cv >= 1 && cv == fix(cv))
                    error('gmm_fit:invalidarg', ...
                        'maxiter should be a positive integer.');
                end
                opts.maxiters = cv;
                
            case 'tol'
                if ~(isfloat(cv) && isreal(cv) && isscalar(cv) && cv > 0)
                    error('gmm_fit:invalidarg', ...
                        'tol should be a real positive scalar.');
                end
                opts.tol = cv;
                
            case 'pricount'
                if ~(isfloat(cv) && isreal(cv) && isscalar(cv) && cv >= 0)
                    error('gmm_fit:invalidarg', ...
                        'pricount should be a real non-negative scalar.');
                end
                opts.pricount = cv;
                
            case 'tied_cov'
                if ~(islogical(cv) && isscalar(cv))
                    error('gmm_fit:invalidarg', ...
                        'tied_cov should be a logical scalar.');
                end
                opts.tied_cov = cv;
                
            case 'display'
                if ~ischar(cv)
                    error('gmm_fit:invalidarg', ...
                        'display should be a string.');
                end
                opts.display = cv;
                
            otherwise
                error('gmm_fit:invalidarg', 'Invalid option name %s', cn);            
        end
    end
end

if isempty(opts.initL)
    opts.initL = randi(K, 1, n);
end





