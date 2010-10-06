function [x, info] = newtonfmin(f, x0, options)
% Perform unconstrained minimization using Newton's method
%
%   x = newtonfmin(f, x0, options);
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
%       x0 is the initial guess of the solution.
%
%       options is the optimization option struct, with following fields:
%
%       - MaxIter:  the maximum number of iterations
%       - TolFun:   the termination tolerance of objective value change
%       - TolX:     the termination tolerance of solution change 
%       - Monitor:  the monitor that shows the procedural information
%
%       The function returns the optimized solution.
%
%   [x, info] = newtonfmin(f, x0, options);
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
%       - Created by Dahua Lin, on Aug 3, 2010
%       - Modified by Dahua Lin, on Aug 7, 2010
%           - use monitor to replace display level.
%

%% verify input 

verify_optim_options('newtonfmin', options, 'MaxIter', 'TolFun', 'TolX');
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

while ~converged && it < options.MaxIter
    
    it = it + 1;
    if omon_level >= optim_mon.IterLevel
        omon.on_iter_start(it);
    end
   
    [v0, g, H] = f(x);
    step = - (H \ g);
    
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



