function [x, info] = gdfmin(f, x0, varargin)
% Perform unconstrained minimization using simple gradient descent
%
%   x = gdfmin(f, x0, ...);
%       solves a (local) minima of f using standard gradient descent.
%       Here, f is the objective function that supports the following 
%       usage:
%
%           [v, g] = f(x);
%
%       In the output, v is the function value evaluated at x, g and H
%       are respectively the gradient and Hessian.
%
%       x0 is the initial guess of the solution.
%
%       One can also specify following options through name/value pairs.
%
%       - MaxIter:  the maximum number of iterations {200}
%       - TolFun:   the termination tolerance of objective value change {1e-6}
%       - TolX:     the termination tolerance of solution change {1e-6}
%       - Display:  the level of display {'none'}|'proc'|'iter'
%
%       The function returns the optimized solution.
%
%   [x, info] = gdfmin(f, x0, ...);
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
    options = struct('MaxIter', 200, 'TolFun', 1e-6, 'TolX', 1e-6);

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


%% main

x = x0;
converged = false;
it = 0;

beta = 0.8;
minstep = 0.5 * options.TolFun;

if omon_level >= optim_mon.ProcLevel
    omon.on_proc_start();
end

total_fcnt = 0;

while ~converged && it < options.MaxIter
    
    it = it + 1;
    if omon_level >= optim_mon.IterLevel
        omon.on_iter_start(it);
    end
    
    [v0, g] = f(x);
    step = -g;
    
    [x, v, dx, fcnt] = linesearch(f, x, v0, step, beta, minstep);
    
    ch = v - v0;
    nrm_dx = norm(dx);
    converged = abs(ch) < options.TolFun && nrm_dx < options.TolX;  
        
    total_fcnt = total_fcnt + (fcnt + 1);
    
    if omon_level >= optim_mon.IterLevel        
        itstat = struct( ...
            'FunValue', v, ...
            'FunChange', ch, ...
            'FunEvals', total_fcnt, ...
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

