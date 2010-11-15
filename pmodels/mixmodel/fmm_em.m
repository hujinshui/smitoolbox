function [FM, Q, status] = fmm_em(gm, X, init, varargin)
% Estimate finite mixture model using E-M algorithm
%
%   [FM, Q] = fmm_em(gm, X, K, ...);
%   [FM, Q] = fmm_em(gm, X, init, ...);
%
%       estimate finite mixture model using E-M algorithm.
%
%       Input arguments:
%       - gm:   the underlying parametric generative model
%       - X:    the observed samples to be fitted
%       - K:    the number of component models
%       - init: the initial guess of component parameters
%       
%       Output arguments:
%       - FM:   the estimated finite mixture model
%       - Q:    the probabilistic assignments of samples [K x n].
%
%       The E-M algorithm iteratively invokes E-steps and M-steps. 
%       In each E-step, the probabilistic assignment of each sample
%       to each component is re-solved; while in E-step, the component
%       parameters are updated based on the new probabilistic assignment.
%
%       If init is given, the function uses it as the initial guess of
%       component parameters, otherwise, K is given, and component
%       parameters are initialized randomly (or using initQ).
%
%       One can further can specify following options in name/value pairs:
%
%       - InitQ:    the initial guess of probabilistic assignment
%                   (default = [], indicating a random initialization)
%
%       - MaxIter:  the maximum number of iterations (default = 100)
%
%       - TolFun:   the tolerance of change of objective function 
%                   at convergence. default = 1e-6
%
%       - TolQ:     the tolerance of change of Q matrix at convergence.
%                   It is measured by L_inf norm. default = 1e-5.
%
%       - Weights:  the sample weights, which can be either a scalar
%                   (shared by all samples) or a vector of length n.
%                   (default = 1).
%
%       - PriCount: the prior counts to each component. 
%                   It can be a scalar, or a vector of length K.
%                   default = 0.
%
%       - Display:  the level of display. It can be either of
%                   'off', 'notify', 'final', 'iter'.
%                   default = 'iter'.
%
%       - InlierConf:   the prior confidence of being inlier.
%                       default = 1.
%
%                       Note that if InlierConf is set to 1. Then
%                       in output, each column of Q sums to 1. 
%                       Otherwise, the sum of each column of Q may
%                       be less than 1, which equals the probability
%                       that the corresponding sample is an inlier.
%
%       - OutlierLogP;  the log-pdf of outlier generative model
%                       default = 0.
%       

%   History
%   -------
%       - Created by Dahua Lin, on Nov 14, 2010
%


%% verify input arguments

if ~isobject(gm)
    error('fmm_em:invalidarg', ...
        'gm should be an object representing a generative model.');
end

n = gm.check_observations(X);
if n <= 0
    error('fmm_em:invalidarg', ...
        'The input observation X is invalid.');
end

if isnumeric(init) && isscalar(init)
    K = init;
    if ~(K == fix(K) && K >= 1)
        error('fmm_em:invalidarg', ...
            'K should be a positive integer.');
    end
    init = [];
else
    K = gm.check_parameters(init);
    if K <= 0
        error('fmm_em:invalidarg', ...
            'The initial guess of parameters is invalid.');
    end
end


% some constants

DISP_NOTIFY = 1;
DISP_FINAL = 2;
DISP_ITER = 3;

opts = struct( ...
    'InitQ', [], ...
    'MaxIter', 100, ...
    'TolFun', 1e-6, ...
    'TolQ', 1e-5, ...
    'Weights', 1, ...
    'PriCount', 0, ...
    'Display', DISP_ITER, ...
    'InlierConf', 1, ...
    'OutlierLogP', 0);

if ~isempty(varargin)
    opts = parlist(opts, varargin{:});
    opts = check_options(opts, K, n);
end

use_robust = opts.InlierConf < 1;


%% main

% initialize status

status = struct( ...
    'niters', 0, ...
    'converged', false, ...
    'objv', NaN, ...    
    'changeF', NaN, ...
    'changeQ', NaN, ...
    'elapsed', 0);

stopTimer = tic;  % start timing

% perform initialization

if isempty(init)
    params = initialize_params(gm, X, n, K, opts.InitQ, opts.Weights);
else
    params = init;
end

if isscalar(opts.PriCount)
    Pi = constmat(K, 1, 1/K);
else
    Pi = (opts.PriCount + 1) / (sum(opts.PriCount) + K);
end

rho = [];
if use_robust    
    l_in = log(opts.InlierConf);
    l_out = log(1 - opts.InlierConf);
    ev_out = constmat(1, n, opts.OutlierLogP) + l_out;
    rho = ones(1, n);
end


ws = opts.Weights;

if opts.Display >= DISP_ITER        
    iter_display(status);    
end

while status.niters < opts.MaxIter
    
    status.niters = status.niters + 1;
    
    if status.niters > 1
        prev_objv = status.objv;
        prev_Q = Q;
    end
    
    % E-step
    
    % evaluate likelihood
    LL = gm.loglik(params, X);    
    
    % infer Q
    Ev = bsxfun(@plus, LL, log(Pi));
    if ~use_robust
        rEv = Ev;
    else
        rEv = bsxfun(@times, rho, Ev);
    end
    Q = nrmexp(rEv, 1);
    
    % infer rho
    if use_robust
        ev_in = dot(Q, Ev, 1) + l_in;
        rho = nrmexp([ev_in; ev_out]);
        rho = rho(1,:);
    end
    
    % Evaluate objective
        
    W = make_w(Q, ws, rho); 
    objv = sum(gm.logpri(params)) ...
        + sum(dot(W, Ev, 1)) ...
        + sum_w(ddentropy(Q), ws);
    
    if ~isequal(opts.PriCount, 0)
        objv = objv + safedot(opts.PriCount, log(Pi));
    end
    if use_robust
        objv = objv ...
            + opts.OutlierLogP * sum_w((1 - rho), ws) ...
            + l_in * sum_w(rho, ws) ...
            + l_out * sum_w(1 - rho, ws) ...
            + sum_w(ddentropy([rho; 1 - rho]), ws);
    end        
    status.objv = objv;
    
    % Determine convergence
    
    if status.niters > 1
        status.changeF = objv - prev_objv;
        status.changeQ = Linfdiff(Q, prev_Q);
        status.converged = ...
            abs(status.changeF) < opts.TolFun && ...
            status.changeQ < opts.TolQ;
        
        if status.converged
            break;
        end
    end       
    
    if opts.Display >= DISP_ITER
        iter_display(status);    
    end
       
    % M-step
    
    % estimate params         
    params = gm.estimate_map(X, W);
    
    % estimate Pi
    Pi = sum(W, 2) + opts.PriCount;
    Pi = Pi / sum(Pi);
        
end

status.elapsed = toc(stopTimer);  % end timing

if opts.Display >= DISP_NOTIFY
    if ~status.converged || opts.Display >= DISP_FINAL
        final_display(status);
    end
end

%% Output

FM = fmm(gm, params, Pi.');
if use_robust
    Q = bsxfun(@times, rho, Q);
end
   
    
%% Auxiliary functions


function W = make_w(Q, ws, rho)

if isempty(ws)
    wm = rho;
else
    if isempty(rho)
        wm = ws;
    else
        wm = ws .* rho;
    end
end

if isempty(wm)
    W = Q;
else
    W = bsxfun(@times, Q, wm);
end



function s = sum_w(x, w)

if isscalar(w)
    s = w * sum(x);
else
    s = dot(x, w);
end



function params = initialize_params(gm, X, n, K, initQ, ws)
    
if isempty(initQ)
    initQ = rand(K, n);
    initQ = bsxfun(@times, initQ, 1 ./ sum(initQ, 1));
end

if isequal(ws, 1)
    wQ = initQ;
elseif isscalar(ws)
    wQ = ws * initQ;
else
    wQ = bsxfun(@times, initQ, ws);
end

params = gm.estimate_map(X, wQ);




function opts = check_options(opts, K, n)

if ~isempty(opts.InitQ)
    q0 = opts.InitQ;
    if ~(isfloat(q0) && ndims(q0) == 2 && size(q0, 1) == K)
        error('fmm_em:invalidopt', ...
            'InitQ should be a numeric matrix with K rows.');
    end
end

mi = opts.MaxIter;
if ~(isscalar(mi) && isnumeric(mi) && mi >= 1 && mi == fix(mi))
    error('fmm_em:invalidopt', ...
        'MaxIter should be a positive integer scalar.');
end

tf = opts.TolFun;
if ~(isscalar(tf) && isfloat(tf) && isreal(tf) && tf > 0)
    error('fmm_em:invalidopt', ...
        'TolFun should be a positive real scalar.');
end

tq = opts.TolQ;
if ~(isscalar(tq) && isfloat(tq) && isreal(tq) && tq > 0)
    error('fmm_em:invalidopt', ...
        'TolQ should be a positive real scalar.');
end

ws = opts.Weights;
if ~(isscalar(ws) || (isvector(ws) && length(ws) == n))
    error('fmm_em:invalidopt', ...
        'Weights should be either a scalar or a real vector of length n.');
end
if size(ws, 1) > 1; opts.Weights = ws.'; end

pc = opts.PriCount;
if ~(isscalar(pc) || (isvector(pc) && length(pc) == K) )
    error('fmm_em:invalidopt', ...
        'PriCount should be a vector of length K.');
end
if size(pc, 2) > 1; opts.PriCount = pc.'; end

dp = opts.Display;
if ~(isnumeric(dp) && isscalar(dp))
    if ~(ischar(dp))
        error('fmm_em:invalidopt', 'Display should be a string.');
    end
    [b, dpi] = ismember(dp, {'off', 'notify', 'final', 'iter'});
    if ~b
        error('fmm_em:invalidopt', 'The value of Display is invalid.');
    end
    opts.Display = dpi - 1;
end

inc = opts.InlierConf;
if ~(isscalar(inc) && isfloat(inc) && isreal(inc) && inc > 0 && inc <= 1)
    error('fmm_em:invalidopt', ...
        'InlierConf should be a real scalar in (0, 1].');
end

olp = opts.OutlierLogP;
if ~(isscalar(olp) && isfloat(olp) && isreal(olp))
    error('fmm_em:invalidarg', ...
        'OutlierLogP should be a real scalar.');
end



%% Display functions

function iter_display(s)

if s.niters == 0
    fprintf('# Iter        objv      objv.ch         Q.ch\n');
    fprintf('-------------------------------------------------------\n');
end

if s.niters >= 1
    fprintf('%4d  %12g %12g %12g\n', s.niters, s.objv, s.changeF, s.changeQ);
end
        

function final_display(s)

fprintf('Terminates with %d iterations (%f secs): ', s.niters, s.elapsed);
if s.converged
    fprintf('converged.\n');
else
    fprintf('NOT converged !\n');
end
    
