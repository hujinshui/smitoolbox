function [x, info] = bfgsfmin(f, x0, varargin)
% Perform unconstrained minimization using (BFGS) quasi-Newton method 
%
%   x = bfgsfmin(f, x0, ...);
%       solves a (local) minima of f using quasi-Newton method, with
%       Broyden-Fletcher-Goldfarb-Shanno (BFGS) updates.
%       Here, f is the objective function that supports the following 
%       usage:
%
%           [v, g] = f(x);
%
%       In the output, v is the function value evaluated at x, g is the 
%       gradient of f at x.
%
%       x0 is the initial guess of the solution.
%
%       One can also specify following options through name/value pairs.
%
%       - MaxIter:  the maximum number of iterations {100}
%       - TolFun:   the termination tolerance of objective value change {1e-8}
%       - TolX:     the termination tolerance of solution change {1e-7}
%       - Display:  the level of display {'none'}|'proc'|'iter'
%       - InitInvHess: the initial inversed Hessian matrix 
%                      (if omitted, the default is the identity matrix)
%
%       The function returns the optimized solution.
%
%   [x, info] = bfgsfmin(f, x0, ...);
%       additionally returns the information of the solution. Here, info
%       is a struct with the following fields:
%
%       - FunValue:     the objective function value
%       - LastChange:   the objective value change at last step
%       - LastMove:     the solution change at last step
%       - IsConverged:  whether the procedure converges
%       - NumIters:     the number of elapsed iterations
%

%   History
%   -------
%       - Created by Dahua Lin, on Aug 7, 2010
%       - Modified by Dahua Lin, on Jan 5, 2010
%           - change the way of option setting
%

%% verify input 

if nargin == 3 && isstruct(varargin{1})
    options = varargin{1};
else
    options = struct('MaxIter', 100, 'TolFun', 1e-8, 'TolX', 1e-7);

    if nargin == 1 && strcmpi(f, 'options')
        x = options;
        return;
    end
    options = smi_optimset(options, varargin{:});
end

omon_level = 0;
if isfield(options, 'Monitor')
    omon = options.Monitor;
    omon_level = omon.level;
end

d = length(x0);
if isfield(options, 'InitInvHess') && ~isempty(options.InitInvHess)
    B = options.InitInvHess;
    if ~(size(B,1) == d && size(B,2) == d)
        error('bfgsfmin:invalidopt', ...
            'InitInvHess should be a d x d real matrix.');
    end
else
    B = eye(d, d);
end


%% main

x = x0;
converged = false;
it = 0;

beta = 0.8;
minstep = 0.5 * options.TolFun;

if omon_level >= optim_mon.ProcLevel
    omon.on_proc_start();
end

while ~converged && it < options.MaxIter
    
    it = it + 1;
    if omon_level >= optim_mon.IterLevel
        omon.on_iter_start(it);
    end
    
    if it == 1
        [v, g] = f(x);
    end
    
    v0 = v;
    g0 = g;
    
    step = - (B * g);    
    [x, v, dx] = linesearch(f, x, v0, step, beta, minstep);
        
    ch = v - v0;
    nrm_dx = norm(dx);
    converged = abs(ch) < options.TolFun && nrm_dx < options.TolX;  
    
    if ~converged  % update B (inverse Hessian)
        [v, g] = f(x);
        y = g - g0;
        sv = y' * dx;
        By = B * y;
        B = B + ((sv + y' * By) / (sv^2)) * (dx * dx') - (1/sv) * (By * dx' + dx * By');
        B = (B + B') / 2;
    end    
    
    if omon_level >= optim_mon.IterLevel        
        itstat = struct( ...
            'FunValue', v, ...
            'FunChange', ch, ...
            'Move', dx, ...
            'MoveNorm', nrm_dx, ...
            'IsConverged', converged);                    
        omon.on_iter_end(it, itstat);
    end
end

if  nargout >= 2 || omon_level >= optim_mon.ProcLevel
    info = struct( ...
        'FunValue', v, ...
        'LastChange', ch, ...
        'LastMove', dx, ...
        'IsConverged', converged, ...
        'NumIters', it);
end

if omon_level >= optim_mon.ProcLevel
    omon.on_proc_end(info);
end

