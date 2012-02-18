function [W, H, objv, converged] = nnmf_std(X, r, W0, varargin)
%NNMF_STD Non-Negative Matrix Factorization
%
%   [W, H] = NNMF_STD(X, r, W0, ...);
%   [W, H, objv, converged] = NNMF_STD( ... );
%
%       Perform Non-negative matrix factorization using standard methods.
%
%       Let X be an m x n matrix, then the goal is to seek two factor
%       matrices W (of size m x r) and H (of size r x n), such that
%       W * H is close to X, under some criterion.
%
%       Input arguments:
%       - X:        The input matrix of size m x n
%       - r:        The rank of the approximation (r <= min(m, n))
%       - W0:       The initial guess of W (empty or an m x r matrix)
%
%       If W0 is left empty, it will be initialized randomly.
%
%       Output arguments:
%       - W:            The solved W-matrix
%       - H:            The solved H-matrix
%       - objv:         The objective value for the solution
%       - converged:    whether the optimization converges
%
%
%       One can specify other options to control the algorithm, in form
%       of name/value pairs:
%
%       - 'method':     The method used to solve the problem
%                       - 'als':    alternate convex optimization
%                                   (only for the case with obj='euc')
%                       - 'mult':   multiplicative update 
%                       default: 'als' when obj is 'euc', 'mult' otherwise.
%       
%       - 'obj':        The objective type
%                       - 'euc':    Euclidean mean squared error (default)
%                       - 'kl':     (average) K-L divergence
%
%       - 'maxiter':    The maximum number of iterations (default = 100)
%
%       - 'tol':        The tolerance of change of W upon convergence.
%                       (default = 1e-6)
%
%       - 'nrm':        The way to normalize the results
%                       - 'WL1':    The L1-norm of each column of W is
%                                   normalized to one (default)
%                       - 'WL2':    The L2-norm of each column of W is
%                                   normalized to one
%
%       - 'display'     The level of display
%                       - 'off':    No display (default)
%                       - 'final':  Only display upon completion
%                       - 'iter':   Display information for each iteration
%

% Created by Dahua Lin, on Feb 17, 2012
%

%% verify input arguments

if ~(isfloat(X) && isreal(X) && ismatrix(X))
    error('nnmf_std:invalidarg', 'X should be a real matrix.');
end
[m, n] = size(X);

if ~(isnumeric(r) && isscalar(r) && r >= 1 && r <= min(m, n))
    error('nnmf_std:invalidarg', ...
        'r should be a positive integer with r <= min(m, n).');
end

if nargin < 3 || isempty(W0)
    W0 = [];
else
    if ~(isfloat(W0) && isreal(W0) && isequal(size(W0), [m r]))
        error('nnmf_std:invalidarg', ...
            'W0 should be a real matrix of size m x r.');
    end
end

[method, obj, maxiter, tol, nrm, displevel] = parse_options(varargin);


%% main

% initialize

if isempty(W0)
    W0 = rand(m, r);
end
wnrms = col_norms(W0, nrm);

W = bsxfun(@times, W0, 1 ./ wnrms);

if method == 2
    H = (W' * W) \ (W' * X);        
    H(H < 0) = 0;
end

% main loop

converged = false;
t = 0;

if displevel < 2
    objv = [];
else
    objv = nan;
end

if displevel >= 2
    fprintf('%6s    %16s (changes)   %8s\n', ...
            'Iter', 'objective', 'ch.W');
    fprintf('--------------------------------------------------------\n');
end

epv = 1e-12;

while ~converged && t < maxiter
    
    t = t + 1;
    Wpre = W;
    
    % do updating
    
    if method == 1      % als
        if obj == 1     % mse            
            H = (W' * W) \ (W' * X);
            H(H < 0) = 0;
            W = (X * H') / (H * H');
            W(W < 0) = 0;
            
        else            % kl
            W = solve_W_kl(X, H);
            H = solve_H_kl(X, W);
            
        end
    else                % mult
        if obj == 1            
            H = H .* ((W' * X) ./ ((W' * W) * H + epv));
            W = W .* (X * H') ./ (W * (H * H') + epv);
            
        else
            C = W' * (X ./ (W * H + epv));
            H = H .* bsxfun(@times, C, 1 ./ sum(W, 1)');
            C = (X ./ (W * H + epv)) * H';
            W = W .* bsxfun(@times, C, 1 ./ sum(H, 2)');
        end            
    end
    
    % normalize
    
    wnrms = col_norms(W, nrm);
    W = bsxfun(@times, W, 1 ./ wnrms);
    H = bsxfun(@times, H, wnrms.');    

    % determine convergence
            
    chW = norm(W - Wpre, inf);
    converged = abs(chW) < tol;          
    
    if displevel >= 2        
        objv_pre = objv;
        objv = eval_objv(X, W, H, obj);
        ch = objv - objv_pre;
                
        fprintf('%6d    %15.6f (%10.3e)   %10.3e\n', ...
            t, objv, ch, chW);
    end
                
end

if nargout >= 3 && isempty(objv)
    objv = eval_objv(X, W, H, obj);
end

if displevel >= 1
    if converged
        fprintf('NNMF converged (with %d iters)\n\n', t);
    else
        fprintf('NNMF did not converge (with %d iters)\n\n', t);
    end
end



%% Auxiliary functions

function v = col_norms(W, nrm)
% calculate per-column norms

if nrm == 1
    v = sum(W, 1);
else
    v = sqrt(sum(W.^2, 1));
end


function v = eval_objv(X, W, H, obj)
% evaluate objective function

Y = W * H;

if obj == 1
    v = norm(X - Y, 'fro')^2 / numel(X);
else
    si = find(X > 0);
    v = sum(X(si) .* (log(X(si)) - log(Y(si)))) + sum(Y(:) - X(:));
    v = v / numel(X);
end


%% Option parsing

function [method, obj, maxiter, tol, nrm, displevel] = parse_options(params)
% parse options

method = [];
obj = 1;
maxiter = 100;
tol = 1e-6;
nrm = 1;
displevel = 0;

if ~isempty(params)
    
    names = params(1:2:end);
    vals = params(2:2:end);
    
    if ~(numel(names) == numel(vals) && iscellstr(names))
        error('nnmf_std:invalidarg', 'Syntax errors for the name/value pairs.');
    end
    
    for i = 1 : numel(names)
        
        cn = names{i};
        cv = vals{i};
        
        switch lower(cn)
            case 'method'
                if ~ischar(cn)
                    error('nnmf_std:invalidarg', ...
                        'The method should be a string');
                end
                switch cv
                    case 'als'
                        method = 1;
                    case 'mult'
                        method = 2;
                    otherwise
                        error('nnmf_std:invalidarg', ...
                            'The method %s is invalid', cv);
                end
                
            case 'obj'
                if ~ischar(cn)
                    error('nnmf_std:invalidarg', ...
                        'The value of obj should be a string');
                end
                switch cv
                    case 'euc'
                        obj = 1;
                    case 'kl'
                        obj = 2;
                    otherwise
                        error('nnmf_std:invalidarg', ...
                            'The value of obj %s is invalid', cv);
                end                                
                
            case 'maxiter'
                if ~(isnumeric(cv) && isscalar(cv) && cv >= 1)
                    error('nnmf_std:invalidarg', ...
                        'The maxiter should be a positive scalar.');
                end
                maxiter = cv;
                
            case 'tol'                                
                if ~(isfloat(cv) && isreal(cv) && isscalar(cv) && cv > 0)
                    error('nnmf_std:invalidarg', ...
                        'The tol should be a positive real scalar.');
                end
                tol = cv;
                
            case 'nrm'                                
                if ~ischar(cn)
                    error('nnmf_std:invalidarg', ...
                        'The value of nrm should be a string');
                end
                switch upper(cv)
                    case 'WL1'
                        nrm = 1;
                    case 'WL2'
                        nrm = 2;
                    otherwise
                        error('nnmf_std:invalidarg', ...
                            'The value of nrm %s is invalid', cv);
                end  
                
            case 'display'
                if ~ischar(cn)
                    error('nnmf_std:invalidarg', ...
                        'The value of display should be a string');
                end
                switch cv
                    case 'off'
                        displevel = 0;
                    case 'final'
                        displevel = 1;
                    case 'iter'
                        displevel = 2;
                    otherwise
                        error('nnmf_std:invalidarg', ...
                            'The value of display: %s is invalid', cv);
                end                  
                
        end                
    end    
    
end


if isempty(method)
    method = obj;
else
    if obj == 2 && method == 1
        error('nnmf_std:invalidarg', ...
            'The method als cannot be used when obj is ''kl''.');
    end
end





