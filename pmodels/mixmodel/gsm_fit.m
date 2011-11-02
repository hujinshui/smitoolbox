function [p, s, converged] = gsm_fit(x, w, K, varargin)
% Fits a Gaussian scale mixture to data
%
%   The pdf of a Gaussian scale mixture (GSM) model is formulated as
%
%       f(x) = \sum_{k=1}^K p_k/sqrt(2 *pi * s_k^2) * 
%                           exp(- x^2 / (2 * s_k^2) );
%
%   Here, p_1, ..., p_K are the prior weights of the components, and
%   s_1, ..., s_K are the scales of the components. Both are parameters
%   to be estimated.
%
%   [p, s] = gsm_fit(x, [], K, ...);
%   [p, s] = gsm_fit(x, w, K, ...);
%       
%       Estimates the parameters of a Gaussian scale mixture, by fitting
%       to the data.
%
%       Input arguments:
%       - x:        a vector of observed values.
%       - w:        the weights of the values. If omitted, each sample
%                   has a weight 1.
%       - K:        the number of mixture components.
%
%       One can specify other options to control the estimation, via
%       name/value pairs:
%
%       - p0:       the initial guess of the component weights p
%
%       - s0:       the initial guess of the component scales s
%
%       - method:   which kind of method to use (default = 'ip')
%                   'ip':   interior-point algorithm
%                   'em':   expectation-maximization algorithm
%
%       - maxiter:  the maximum number of iterations
%                   default = 100
%
%       - iterlen:  the number of E-M steps in an iteration 
%                   (default = 5)
%
%       - tol:      the tolerance of change of objective at convergence
%                   (default = 1e-6)
%
%       - tolx:     the tolerance of change of solution at convergence
%                   (default = 1e-6)
%
%       - display:  the level of display, which can be 'none', 'notify',
%                   'final', and 'iter'.
%
%       If p0 and s0 are not given, then this function uses its own
%       method to give an initial guess.
%
%   [p, s, converged] = gsm_fit( ... );
%
%       returns whether the optimization procedure converges within
%       tolerance.
%
%

% Created by Dahua Lin, on Nov 2, 2011
%

%% verify input arguments

if ~(isfloat(x) && isvector(x) && isreal(x) && ~issparse(x))
    error('gsm_fit:invalidarg', 'x should be a non-sparse real vector.');
end
n = numel(x);
if size(x, 1) > 1; x = x .'; end  % turn x into a row vector

if ~isempty(w)
    if ~(isfloat(w) && isvector(w) && isreal(w) && ~issparse(w))
        error('gsm_fit:invalidarg', 'w should be a non-sparse real vector.');
    end
    if numel(w) ~= n
        error('gsm_fit:invalidarg', 'The size of w is not consistent with x.');
    end
    if size(w, 2) > 1; w = w.'; end % turn w into a column vector
else
    w = ones(n, 1);
end

if ~(isscalar(K) && isnumeric(K) && K == fix(K) && K >= 1)
    error('gsm_fit:invalidarg', 'K should be a positive integer number.');
end

% check options

[p0,s0,method,maxiter,iterlen,tol,tolx,dispstr] = check_options(varargin);


%% main skeleton

% initialization

x2 = x.^2;

if isempty(p0)
    p0 = constmat(K, 1, 1.0 / double(K));
else
    p0 = p0 / sum(p0); % ensure that p0 sums to one
end

if isempty(s0)
    beta0 = init_beta0(x2, p0);
else
    beta0 = 1 ./ (s0.^2);
end

switch method
    case 'ip'
        [p, beta, converged] = fit_ip(x2, w, p0, beta0, ...
            maxiter, tol, tolx, dispstr);
    case 'em'
        [p, beta, converged] = fit_em(x2, w, p0, beta0, ...
            maxiter, iterlen, tol, tolx, dispstr);
end

s = sqrt(1 ./ beta);


%% numeric solver 

function [p, beta, converged] = fit_ip(x2, w, p0, beta0, maxiter, tol, tolx, dispstr)

K = numel(beta0);
fun = @(y) direct_objfun(y, x2, w, K);

sol0 = [p0; beta0];

Aeq = [ones(1, K), zeros(1, K)];
lb = zeros(2*K, 1);

options = optimset('Algorithm', 'interior-point', ...
    'MaxIter', maxiter, ...
    'TolFun', tol, ...
    'TolX', tolx, ...
    'Display', dispstr);

[sol, ~, eflag] = fmincon(fun, sol0, [], [], Aeq, 1, lb, [], [], options);

p = sol(1:K);
beta = sol(K+1:end);
converged = (eflag > 0);


% objective function

function [v, g] = direct_objfun(sol, x2, w, K)

p = sol(1:K);
beta = sol(K+1:end);

U = (-0.5 * beta) * x2;  % outer product

lrho = log(p) + 0.5 * log(beta);

E = bsxfun(@plus, lrho, U);

maxE = max(E, [], 1);
rE = bsxfun(@minus, E, maxE);
rE = exp(rE);
rE_sum = sum(rE, 1);

Lik = log(rE_sum) + maxE;
v = Lik * w;
v = -v;  % we are to maximize

if nargout >= 2     % evaluate gradient
    
    Q = bsxfun(@times, rE, 1 ./ rE_sum);
    
    gp = Q * w;
    gbeta = 0.5 * (gp ./ beta - Q * (w .* (x2.')));

    g = [gp; gbeta];
    g = -g; 
end


%% expectation maximization

function [p, beta, converged] = fit_em(x2, w, p0, beta0, maxiter, iterlen, tol, tolx, dispstr)

displevel = get_displevel(dispstr);

p = p0;
beta = beta0;

wx2 = w .* x2.';

it = 0;
converged = false;

if displevel >= 3
    fprintf('%6s    %12s    %12s    %12s\n', 'Iter', 'objective', 'objv-change', 'sol-change');
end


E = bsxfun(@minus, log(p) + log(beta) * 0.5, (0.5 * beta) * x2);

while ~converged && it < maxiter

    it = it + 1;
    
    for t = 1 : iterlen
    
        % E-step
        
        Q = nrmexp(E, 1);
        
        % M-step
        
        if it > 1
            pre_p = p;
            pre_beta = beta;
        end
        
        sw = Q * w;
        beta = sw ./ (Q * wx2);
        p = sw / sum(sw);    
        
        E = bsxfun(@minus, log(p) + log(beta) * 0.5, (0.5 * beta) * x2);
    end
    
    % Evaluate objective
        
    ent = ddentropy(Q);
    
    if it > 1
        prev_objv = objv;
    end
    objv = (sum(Q .* E, 1) + ent) * w;
        
    if it == 1
        ch_fun = nan;
        ch_sol = nan;
    else        
        ch_fun = objv - prev_objv;
        ch_sol = max( norm(p - pre_p, inf), norm(beta - pre_beta, inf) );
    end
    
    converged = it > 1 && abs(ch_fun) < tol && ch_sol < tolx;
    
    % Display
    
    if displevel >= 3
        fprintf('%6d    %12.5g    %12.5g    %12.5g\n', it, objv, ch_fun, ch_sol);
    end
end

if displevel >= 2 || (displevel >= 1 && ~converged)
    if converged
        fprintf('GSM fitting (by EM) converged (with %d iterations)\n', it);
    else
        fprintf('GSM fitting (by EM) have not converged (with %d iterations)\n', it);
    end
end
  

%% auxiliary functions

function beta0 = init_beta0(x2, p0)

n = numel(x2);
sx2 = sort(x2);

ei = round(cumsum(p0) * n);
ei(end) = n;
si = [1; ei(1:end-1) + 1];

K = numel(p0);
beta0 = zeros(K, 1);
for k = 1 : K
    cx2 = sx2(si(k):ei(k));
    cn = ei(k) - si(k) + 1;
    beta0(k) = cn / sum(cx2);
end


function [p0,s0,method,maxiter,iterlen,tol,tolx,dispstr] = check_options(nvlist)

p0 = [];
s0 = [];
method = 'ip';
maxiter = 100;
iterlen = 5;
tol = 1e-6;
tolx = 1e-6;
dispstr = 'none';

if ~isempty(nvlist)
    
    onames = nvlist(1:2:end);
    ovals = nvlist(2:2:end);
    
    if ~(numel(onames) == numel(ovals) && iscellstr(onames))
        error('gsm_fit:invalidarg', 'The name/value list is invalid.');
    end
    
    for i = 1 : numel(onames)        
        cn = onames{i};
        cv = ovals{i};
        
        switch lower(cn)
            case 'p0'
                if ~(isfloat(cv) && isreal(cv) && isequal(size(cv), [K 1]))
                    error('gsm_fit:invalidarg', ...
                        'p0 should be a K x 1 real valued vector.');
                end
                p0 = cv;
                
            case 's0'                                
                if ~(isfloat(cv) && isreal(cv) && isequal(size(cv), [K 1]))
                    error('gsm_fit:invalidarg', ...
                        's0 should be a K x 1 real valued vector.');
                end
                s0 = cv;
                
            case 'method'
                if ~(ischar(cv) && (strcmpi(cv, 'ip') || strcmpi(cv, 'em')))
                    error('gsm_fit:invalidarg', ...
                        'The method must be a string, which can be ''ip'' or ''em''.');
                end
                method = lower(cv);
                
            case 'maxiter'
                if ~(isnumeric(cv) && isscalar(cv) && cv >= 1)
                    error('gsm_fit:invalidarg', ...
                        'maxiter must be a positive number.');
                end
                maxiter = cv;
                                            
            case 'iterlen'
                if ~(isnumeric(cv) && isscalar(cv) && cv >= 1)
                    error('gsm_fit:invalidarg', ...
                        'iterlen must be a positive number.');
                end
                iterlen = cv;                
                
            case 'tol'
                if ~(isfloat(cv) && isscalar(cv) && isreal(cv) && cv > 0)
                    error('gsm_fit:invalidarg', ...
                        'tol must be a positive real value.');
                end
                tol = cv;
                                            
            case 'tolx'
                if ~(isfloat(cv) && isscalar(cv) && isreal(cv) && cv > 0)
                    error('gsm_fit:invalidarg', ...
                        'tolx must be a positive real value.');
                end
                tolx = cv;                
                
            case 'display'
                if ~(ischar(cv) && ismember(cv, {'none', 'notify', 'final', 'iter'}))
                    error('gsm_fit:invalidarg', ...
                        'The value of display is invalid.');
                end
                dispstr = cv;
                
        end
    end
    
end



function displevel = get_displevel(dispstr)

switch dispstr
    case 'none'
        displevel = 0;
    case 'notify'
        displevel = 1;
    case 'final'
        displevel = 2;
    case 'iter'
        displevel = 3;
end






