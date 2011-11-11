function [x, info] = newtonfmin(f, x0, varargin)
% Perform unconstrained minimization using Newton's method
%
%   x = newtonfmin(f, x0, ...);
%       solves a (local) minima of f using standard Newton's method.
%       Here, f is the objective function that supports the following 
%       usage:
%
%           [v, g, H] = f(x);
%
%       In the output, v is the function value evaluated at x, g and H
%       are respectively the gradient and Hessian.
%
%       Note that H need not to be a true Hessian matrix, the function
%       just uses it as H \ g to derive the updating direction. Therefore,
%       any object that supports left matrix division is fine.
%
%       If the option 'DirectNewton' is set to true, then f should support
%
%           [v, dx] = f(x);
%
%       Here, dx is the Newton direction (-H\g) evaluated at x. 
%       This is useful when there is some structure of H and g such that
%       dx can be evaluated more efficiently with special implementation.
%
%       x0 is the initial guess of the solution.
%
%       One can also specify following options through name/value pairs.
%
%       - MaxIter:  the maximum number of iterations {30}
%       - TolFun:   the termination tolerance of objective value change {1e-9}
%       - TolX:     the termination tolerance of solution change {1e-7}
%       - Display:  the level of display {'none'}|'proc'|'iter'
%       - DirectNewton: whether f directly outputs the newton direction
%                       (-H\g) as the second argument. {false}
%
%       The function returns the optimized solution.
%
%   [x, info] = newtonfmin(f, x0, ...);
%       additionally returns the information of the solution. Here, info
%       is a struct with the following fields:
%
%       - FunValue:     the objective function value
%       - LastChange:   the objective value change at last step
%       - LastMove:     the solution change at last step
%       - IsConverged:  whether the procedure converges
%       - NumIters:     the number of elapsed iterations
%
%   options = newtonfmin('options');
%       gets the default options struct.
%

%   History
%   -------
%       - Created by Dahua Lin, on Aug 3, 2010
%       - Modified by Dahua Lin, on Aug 7, 2010
%           - use monitor to replace display level.
%       - Modified by Dahua Lin, on Jan 5, 2010
%           - change the way of option setting
%       - Modified by Dahua Lin, on April 23, 2011
%           - support DirectNewton.
%

%% check options

if nargin == 3 && isstruct(varargin{1})
    options = varargin{1};
else
    options = struct('MaxIter', 30, 'TolFun', 1e-9, 'TolX', 1e-7);

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

if isfield(options, 'DirectNewton')
    dnewton = options.DirectNewton;
else
    dnewton = false;
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
   
    if ~dnewton        
        [v0, g, H] = f(x);
        step = - (H \ g);
    else
        [v0, step] = f(x);
    end
    
    [x, v, dx] = linesearch(f, x, v0, step, beta, minstep);
    
    ch = v - v0;
    nrm_dx = norm(dx);
    converged = abs(ch) < options.TolFun && nrm_dx < options.TolX;  
        
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



