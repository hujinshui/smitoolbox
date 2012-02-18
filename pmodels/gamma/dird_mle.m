function [alpha, objv] = dird_mle(X, w, alpha0, varargin)
%DIRD_MLE Maximum-likelihood Estimation of Dirichlet Distribution
%
%   [alpha, objv] = DIRD_MLE(X, w, alpha0);
%   [alpha, objv] = DIRD_MLE(X, w, alpha0, ...);
%
%       Performs maximum-likelihood estimation of Dirichlet distribution.
%
%       Input arguments
%       ---------------
%           X:              The observed samples or their statistics [d x n]
%           w:              The sample weights [empty or 1 x n]
%           alpha0:         The initial guess of alpha [empty or d x 1]
%                           If empty, alpha0 will be initialized to 
%                           1 + rand(d, 1).
%
%           If the 4th argument is 'normal', then X is treated as normal
%           samples drawn from the Dirichlet distribution. 
%
%           If the 4th argument is 'stat', then X is considered to be
%           (the expectation) of log(samples). This syntax can be 
%           useful in various context, e.g. the E-M algorithm for latent 
%           Dirichlet allocation.
%
%           If the 4th argument is not given, it is set to 'normal'.
%           
%       Output arguments
%       ----------------
%           alpha:          The solved alpha vectors
%           objv:           The objective value of the solution.
%
%       One can input other parameters to control the optimization, 
%       in form of name/value pairs
%
%           - input:        'normal' or 'stat':
%                           'normal': X are samples from the Distribution
%                           'stat':   X are expectation of log(samples)
%                           default = 'normal'.
%
%           - maxiter:      The maximum number of iterations {100}
%           - tolfun:       The tolerance of function value changes {1e-9}
%           - tolx:         The tolerance of solution changes {1e-8}
%           - display:      the level of display {'off'}|'notify'|'final'|'iter'
%
%       Note: fmincon will be used to optimize the function, and these
%       options will be input to newtonfmin.
%       

%% parse inputs

if ~(isfloat(X) && isreal(X) && ismatrix(X))
    error('dird_mle:invalidarg', 'X should be a real matrix.');
end
[d, n] = size(X);

if ~isempty(w)
    if ~(isfloat(w) && isreal(w) && isvector(w) && numel(w) == n)
        error('dird_mle:invalidarg', 'w should be a vector of length n.');
    end
end

if isempty(alpha0)
    alpha0 = 1 + rand(d, 1);
else
    if ~(isfloat(alpha0) && isreal(alpha0) && isequal(size(alpha0), [d 1]))
        error('dird_mle:invalidarg', 'alpha0 should be a d x 1 real vector.');
    end    
end

[is_normal, opts] = parse_options(varargin);


%% main

% compute stats

if is_normal
    v = log(X);
else
    v = X;
end
if n > 1
    if isempty(w)
        tw = n;
        v = sum(v, 2) * (1 / n);
    else
        tw = sum(w);
        v = v * (w(:) / tw);
    end
else
    if isempty(w)
        tw = 1;
    else
        tw = w;
    end
end
    

% optimize

objfun = @(a) dird_mle_objfun(a, v);

[alpha, fv] = fmincon(objfun, alpha0, [], [], [], [], ...
    zeros(d, 1), [], [], opts);

objv = (-fv * tw);


%% objective function

function [v, g, H] = dird_mle_objfun(a, t)

sa = sum(a);
v = sum(gammaln(a)) - gammaln(sa) - sum((a - 1) .* t, 1);

if nargout >= 2
    g = psi(a) - psi(sa) - t;
end

if nargout >= 3
    H = diag(psi(1, a)) - psi(1, sa);     
end



%% Function parsing

function [is_normal, opts] = parse_options(params)

is_normal = 1;
maxiter = 100;
tolfun = 1e-9;
tolx = 1e-8;
display = 'off';

if ~isempty(params)
    
    names = params(1:2:end);
    vals = params(2:2:end);
    
    if ~(numel(names) == numel(vals) && iscellstr(names))
        error('dird_mle:invalidarg', 'Invalid name/value list.');
    end
    
    for i = 1 : numel(names)
        
        cn = names{i};
        v = vals{i};
        
        switch lower(cn)
            
            case 'input'
                if strcmp(v, 'normal')
                    is_normal = 1;
                elseif strcmp(v, 'stat')
                    is_normal = 0;
                else
                    error('dird_mle:invalidarg', ...
                        'The value of option input should be ''normal'' or ''stat''.');
                end
        
            case 'maxiter'
                if ~(isnumeric(v) && isscalar(v) && v >= 1)
                    error('dird_mle:invalidarg', ...
                        'maxtter should be a positive integer scalar.');
                end
                maxiter = v;
                
            case 'tolfun'
                if ~(isfloat(v) && isscalar(v) && v > 0)
                    error('dird_mle:invalidarg', ...
                        'tolfun should be a positive scalar.');
                end
                tolfun = v;
                
            case 'tolx'
                if ~(isfloat(v) && isscalar(v) && v > 0)
                    error('dird_mle:invalidarg', ...
                        'tolx should be a positive scalar.');
                end
                tolx = v;
                
            case 'display'
                if ~(ischar(v) || ...
                        ismember(v, {'off', 'notify', 'final', 'iter'}))
                    error('dird_mle:invalidarg', ...
                        'The value of display is invalid.');
                end
                display = v;

            otherwise
                error('dird_mle:invalidarg', ...
                    'The option name %s is unknown.', cn);                
        end    
    end
    
end

opts = optimset( ...
    'MaxIter', maxiter, ...
    'TolFun', tolfun, ...
    'TolX', tolx, ...
    'Display', display, ...
    'GradObj', 'on', ...
    'Hessian', 'user-supplied');






