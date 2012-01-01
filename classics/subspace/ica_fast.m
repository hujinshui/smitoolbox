function W = ica_fast(X, k, W0, varargin)
%ICA_FAST Fast Independent Component Analysis (Fast ICA)
%
%   W = ICA_FAST(X, k);
%   W = ICA_FAST(X, k, [], ...);
%   W = ICA_FAST(X, k, W0, ...);
%
%       Performs Fast ICA on a given set of data. 
%
%       Input arguments:
%       - X:        the sample matrix, where each column is a sample.
%       - k:        the number of components to derive.
%       - W0:       the initial guess of the weight matrix W.
%
%       Output argument:
%       - W:        the matrix comprised of weight vectors. Each column
%                   of W is a weight vector, and the size of W is d x k.
%
%       One can specify other options in form of name/value pairs to
%       control the algorithm.
%
%       - 'method'      The method to solve multiple weight vectors.
%                       - 'defl':   choose the vectors one by one
%                       - 'symm':   perform symmetric decorrelation
%                       (default is 'defl')
%
%       - 'g'           the name of the nonlinearity function g:
%                       - 'pow3':       g(u)=u^3
%                       - 'tanh':       g(u)=tanh(a*u)
%                       - 'gauss':      g(u)=u*exp(-a*u^2/2)
%                       - 'skew':       g(u)=u^2
%                       The default choice is 'pow3'.
%
%       - 'gparam'      the parameter (a) of the nonlinearity function.
%                       (default value = 1)
%
%       - 'tol'         the tolerance of change at convergence.
%                       (default = 1e-6)
%
%       - 'maxiters'    the maximum number of iterations (default = 200)
%
%       - 'verbose'     whether to show procedural information.
%                       (default = false)
%       

% Created by Dahua Lin, on Dec 31, 2012
%

%% verify input arguments

if ~(isfloat(X) && isreal(X) && ndims(X) == 2)
    error('ica_fast:invalidarg', 'X should be a real matrix.');
end
d = size(X, 1);

if ~(isnumeric(k) && isscalar(k) && k == fix(k) && k <= d)
    error('ica_fast:invalidarg', 'k should be a positive integer with k <= d.');
end

if nargin >= 3 && ~isempty(W0)
    if ~(isfloat(W0) && isreal(W0) && isequal(size(W0), [d k]))
        error('ica_fast:invalidarg', 'W0 should be a d x k real matrix.');
    end
else
    W0 = [];
end

[use_symm, gfun, tol, maxiters, verbose] = getopts(varargin);


%% main

% initialize

if isempty(W0)
    W = orth(randn(d, k));
else
    W = W0;
end

% main loop

if verbose
    disp('Fast ICA Running:');
end

if k == 1
    if verbose
        fprintf('Solving single component ...\n');
    end
    W = defl_optimize(X, W, [], gfun, maxiters, tol, verbose);
else
    if use_symm  
        if verbose
            fprintf('Solving all components with symmetric decorrelation ...\n');
        end
        W = symm_optimize(X, W, gfun, maxiters, tol, verbose);        
        
    else   
        if verbose
            fprintf('Solving component %d ...\n', 1);
        end
        W(:,1) = defl_optimize(X, W(:,1), [], gfun, maxiters, tol, verbose);
        for i = 2 : k
            fprintf('Solving component %d ...\n', i);
            W(:,i) = defl_optimize(X, W(:,i), W(:,1:i-1), gfun, ...
                maxiters, tol, verbose);
        end
        
    end
end

if verbose
    disp('Fast ICA finished.');
end


%% Core optimization routines

function w = defl_optimize(X, w, B, gfun, maxiters, tol, verbose)
% optimizing a particular w

n = size(X, 2);
t = 0;
converged = false;

if ~isempty(B)
    w = w - B * (B' * w);    
end
w = w * (1 / norm(w));

while ~converged && t < maxiters
    t = t + 1;
    w_pre = w;
    
    [v, f] = gfun(w' * X);
    w = (X * v') * (1/n) - (sum(f) / n) * w;
    
    if ~isempty(B)
        w = w - B * (B' * w);
    end
    w = w * (1 / norm(w));
    
    dev = norm(w - w_pre);
    converged = dev < tol;
    if verbose
        fprintf('\tIter %d: change = %.4g\n', t, dev);
    end
end

if verbose
    if converged
        fprintf('\tConverged.\n');
    else
        fprintf('\tNOT converged.\n');
    end
end



function W = symm_optimize(X, W, gfun, maxiters, tol, verbose)
% optimizing the whole W in parellel

n = size(X, 2);
t = 0;
converged = false;

W = bsxfun(@times, W, 1 ./ sqrt(sum(W.^2, 1)));

while ~converged && t < maxiters
    t = t + 1;
    W_pre = W;
    
    [V, F] = gfun(W' * X);
    W = (X * V') * (1/n) - bsxfun(@times, W, (sum(F,2) / n).');
    
    [U,S,V] = svd(W, 0); %#ok<ASGLU>
    W = U * V';

    dev = max(sqrt(sum((W - W_pre).^2, 1)));
    converged = dev < tol;
    
    if verbose
        fprintf('\tIter %d: change = %.4g\n', t, dev);
    end
end

if verbose
    if converged
        fprintf('\tConverged.\n');
    else
        fprintf('\tNOT converged.\n');
    end
end



%% nonlinearity functions (g)

function [V, F] = nf_pow3(Y)

V = Y.^3;
F = 3 * (Y.^2);


function [V, F] = nf_tanh(Y, a)

if a == 1
    AY = Y;
else
    AY = a * Y;
end
V = tanh(AY);
F = a ./ (cosh(AY) .^ 2);


function [V, F] = nf_gauss(Y, a)

E = exp((-a/2) * (Y.^2));
V = Y .* E;
F = (1 - a * (Y.^2)) .* E;


function [V, F] = nf_skew(Y)

V = Y.^2;
F = 2 * Y;



%% option parsing

function [use_symm, gfun, tol, maxiters, verbose] = getopts(params)

use_symm = false;
g = 'pow3';
a = 1;
tol = 1e-6;
maxiters = 200;
verbose = false;

if ~isempty(params)
    
    onames = params(1:2:end);
    ovals = params(2:2:end);
    
    for i = 1 : numel(onames)
        cn = onames{i};
        cv = ovals{i};
        
        switch cn
            case 'method'
                if ischar(cv)
                    if strcmpi(cv, 'defl')
                        use_symm = false;
                    elseif strcmpi(cv, 'symm')
                        use_symm = true;
                    else
                        error('ica_fast:invalidarg', ...
                            'Invalid method %s', cv);
                    end
                else
                    error('ica_fast:invalidarg', ...
                        'The method should be a string.');
                end
                
            case 'g'
                if ischar(cv)
                    g = lower(cv);
                else
                    error('ica_fast:invalidarg', ...
                        'The value of the option g should be a string.');
                end
                
            case 'gparam'
                if ~(isfloat(cv) && isscalar(cv) && isreal(cv))
                    error('ica_fast:invalidarg', ...
                        'The value of gparam should be a real scalar.');
                end
                a = cv;
                
            case 'tol'
                if ~(isfloat(cv) && isscalar(cv) && isreal(cv) && cv > 0)
                    error('ica_fast:invalidarg', ...
                        'The value of tol should be a positive real scalar.');
                end
                tol = cv;
                
            case 'maxiters'
                if ~(isnumeric(cv) && isscalar(cv) && cv >= 1)
                    error('ica_fast:invalidarg', ...
                        'The value of maxiters should be a number no less than 1.');
                end
                maxiters = cv;
                
            case 'verbose'
                if ~(islogical(cv) && isscalar(cv))
                    error('ica_fast:invalidarg', ...
                        'The value of verbose should be a logical scalar.');
                end
                verbose = cv;            
        end                
    end    
end

switch g
    case 'pow3'
        gfun = @nf_pow3;
    case 'tanh'
        gfun = @(x) nf_tanh(x, a);
    case 'gauss'
        gfun = @(x) nf_gauss(x, a);
    case 'skew'
        gfun = @(x) nf_skew(x);
    otherwise
        error('ica_fast:invalidarg', ...
            'Invalid g-function name %s', g);
end
        
