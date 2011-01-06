function [a, rw, info] = irls(X, y, wf, s, r2, a0, varargin)
% Iterative Reweighted Least Square
%
%   The iterative reweighted least square algorithm updates vector a
%   by alternating between the following two steps:
%
%       1. solve the weighted least square problem to get a
%
%           minimize sum_i w_i || x_i' * a - y_i ||^2 + (r2/2) * ||a||^2
%
%       2. update the weights as
%
%           w_i = wf(||x_i' * a - y_i|| / s)
%
%          Here, wf is a function that calculates a weight based on 
%          the residue norm.
%
%   a = irls(X, y, wf, s, r2, a0, ...);
%   [a, rw] = irls(X, y, wf, s, r2, a0, ...);
%   [a, rw, info] = irls(X, y, wf, s, r2, a0, ...);
%
%       performs iterative reweighted least square in solving the
%       coefficient vector/matrix a.
%
%       Input arguments:
%       - X:        the design matrix. In particular, x_i is given by
%                   the i-th row of X, i.e. X(i,:).
%
%       - y:        the response vector/matrix. y(i,:) corresponds to 
%                   X(i,:).
%
%       - wf:       the weighting function, which should be able to take
%                   into multiple input residues in form of a matrix, and
%                   return a matrix of the same size.
%
%       - s:        the scale parameter
%
%       - a0:       the initial guess of a.
%
%       Output arguments:
%       - a:        the resultant coefficient vector/matrix.
%
%       - rw:       the weights used in the final iteration (n x 1)
%
%       - info:     a struct that contains the procedural information.
%
%       Suppose X is a matrix of size n x d, and y is a matrix of size
%       n x q, then a (and a0) will be a matrix of size d x q.
%
%       One can specify additional options to control the procedure in
%       name/value pairs.
%
%       - MaxIter:      the maximum number of iterations {100}
%       - TolFun:       the tolerance of objective value change at
%                       convergence {1e-8}
%       - TolX:         the tolerance of change of a at convergence {1e-8}
%       - Display:      the level of information displaying
%                       {'none'|'proc'|'iter'}
%       - Monitor:      the monitor that responses to procedural updates
%

%   History
%   -------
%       - Created by Dahua Lin, on Jan 6, 2011
%

%% verify input and check options

if ~(isfloat(X) && ndims(X) == 2)
    error('irls:invalidarg', 'X should be a numeric matrix.');
end
[n, d] = size(X);

if ~(isfloat(y) && ndims(y) == 2 && size(y, 1) == n)
    error('irls:invalidarg', 'y should be a numeric matrix with n rows.');
end
q = size(y, 2);

if ~isa(wf, 'function_handle')
    error('irls:invalidarg', 'wf should be a function handle.');
end

if ~(isfloat(s) && isscalar(s) && s > 0)
    error('irls:invalidarg', 's should be a positive scalar.');
end

if ~(isfloat(r2) && isscalar(r2) && r2 >= 0)
    error('irls:invalidarg', 'r2 should be a non-negative scalar.');
end

if ~(isfloat(a0) && isequal(size(a0), [d q]))
    error('irls:invalidarg', 'a0 should be a numeric matrix of size d x q.');
end

if numel(varargin) == 1 && isstruct(varargin{1})
    options = varargin{1};
else
    options = struct('MaxIter', 100, 'TolFun', 1e-8, 'TolX', 1e-8);

    if nargin == 1 && strcmpi(f, 'options')
        a = options;
        return;
    end
    options = smi_optimset(options, varargin{:});
end

omon_level = 0;
if isfield(options, 'Monitor')
    omon = options.Monitor;
    omon_level = omon.level;
end


%% main

a = a0;
converged = false;
it = 0;

if omon_level >= optim_mon.ProcLevel
    omon.on_proc_start();
end

% initial weighting

e = X * a - y;
[rn, rn2] = e_to_rn(e);
rw = wf(rn);
v = (rw' * rn2) / 2;


while ~converged && it < options.MaxIter
    
    it = it + 1;
    if omon_level >= optim_mon.IterLevel
        omon.on_iter_start(it);
    end
   
    a_p = a;
    v_p = v;
    
    % re-solve a
    a = llsq(X, y, rw, r2);
    
    % re-evaluate rw and v
    e = (1/s) * (X * a - y);
    [rn, rn2] = e_to_rn(e);
    rw = wf(rn);
    v = (rw' * rn2) / 2;
        
    % determine convergence
    ch = v - v_p;
    nrm_da = norm(a - a_p);
    converged = abs(ch) < options.TolFun && nrm_da < options.TolX;  
        
    
    if omon_level >= optim_mon.IterLevel        
        itstat = struct( ...
            'FunValue', v, ...
            'FunChange', ch, ...
            'Move', NaN, ...
            'MoveNorm', nrm_da, ...
            'IsConverged', converged);                    
        omon.on_iter_end(it, itstat);
    end
end


if  nargout >= 2 || omon_level >= optim_mon.ProcLevel
    info = struct( ...
        'FunValue', v, ...
        'LastChange', ch, ...
        'LastMove', nrm_da, ...
        'IsConverged', converged, ...
        'NumIters', it);
end

if omon_level >= optim_mon.ProcLevel
    omon.on_proc_end(info);
end


function [rn, rn2] = e_to_rn(e)

if size(e, 2) == 1
    rn = abs(e);
    rn2 = e .^ 2;
else
    rn2 = dot(e, e, 2);
    rn = sqrt(rn2);
end


