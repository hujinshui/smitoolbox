function S = fmm_init(pri, gm, X, w, method, arg)
%FMM_INIT Initialize an FMM solution
%
%   S = FMM_INIT(pri, gm, X, w, 'params', A);
%
%       Initializes a finite mixture model solution with initial
%       parameters (given by A).
%
%       Other inputs:
%       - pri:      The prior object
%       - gm:       The generative model
%       - X:        The observations
%       - w:        The sample weights
%
%   S = FMM_INIT(pri, gm, X, w, 'labels',  z);
%
%       Initializes a finite mixture model with initial labels.
%
%   S = FMM_INIT(pri, gm, X, w, 'Q', Q);
%
%       Initializes a finite mixture model with the initial soft
%       assignment matrix.
%
%   S = FMM_INIT(pri, gm, X, w, 'rand', K);
%
%       Randomly initializes a finite mixture model with K components.
%

% Created by Dahua Lin, on Dec 27, 2011
%

%% verify arguments

if ~(isfloat(X) && isreal(X) && ndims(X) == 2)
    error('fmm_init:invalidarg', ...
        'The observation matrix X should be a real matrix.');
end
n = size(X, 2);

if isempty(w)
    w = [];
else
    if ~(isfloat(w) && isreal(w) && isvector(w) && numel(w) == n)
        error('fmm_init:invalidarg', ...
            'w should be a real vector of length n.');
    end
    if size(w, 2) > 1
        w = w.';
    end
end

if ~ischar(method)
    error('fmm_init:invalidarg', 'The method must be a string.');
end

%% main delegate

switch lower(method)
    case 'params'
        [K, params] = fmm_init_params(pri, gm, X, w, arg);
        
    case 'labels'
        [K, params] = fmm_init_labels(pri, gm, X, w, arg);
        
    case 'q'
        [K, params] = fmm_init_Q(pri, gm, X, w, arg);
        
    case 'rand'
        [K, params] = fmm_init_rand(pri, gm, X, w, arg);
        
    otherwise
        error('fmm_init:invalidarg', ...
            'Unknown method name %s', method);
end

S.K = K;
S.Pi = constmat(K, 1, 1.0/K);
S.params = params;


%% core functions

function [K, params] = fmm_init_params(pri, gm, X, w, params) %#ok<INUSL>

K = gm.query_params(params);


function [K, params] = fmm_init_labels(pri, gm, X, w, z)

if ~(isnumeric(z) && ~issparse(z) && isvector(z) && isreal(z))
    error('fmm_init:invalidarg', 'z should be a real vector.');
end
if size(z, 1) > 1; z = z.'; end

K = max(z);
params = fmm_est_params(pri, gm, X, w, {K, z});


function [K, params] = fmm_init_Q(pri, gm, X, w, Q)

if ~(isfloat(Q) && isreal(Q) && ndims(Q) == 2)
    error('fmm_init:invalidarg', 'Q should be a real matrix.');
end
K = size(Q, 1);
params = fmm_est_params(pri, gm, X, w, Q);


function [K, params] = fmm_init_rand(pri, gm, X, w, K)

n = gm.query_obs(X);
if ~(isnumeric(K) && isscalar(K) && K == fix(K) && K >= 1 && K <= n)
    error('fmm_init:invalidarg', 'K should be a positive integer with K <= n.');
end

z = ceil((K/n) * (1:n));
z(end) = K;

% random suffle
[~, si] = sort(rand(1, n));
z = z(si);

params = fmm_est_params(pri, gm, X, w, {K, z});

