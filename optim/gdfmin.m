function [x, info] = gdfmin(f, x0, options)
% Perform unconstrained minimization using simple gradient descent
%
%   x = gdfmin(f, x0, options);
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
%       options is the optimization option struct, with following fields:
%
%       - MaxIter:  the maximum number of iterations
%       - TolFun:   the termination tolerance of objective value change
%       - TolX:     the termination tolerance of solution change 
%       - Monitor:  the monitor that shows the procedural information
%
%       The function returns the optimized solution.
%
%   [x, info] = gdfmin(f, x0, options);
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
%

%% verify input 

verify_optim_options('gdfmin', options, 'MaxIter', 'TolFun', 'TolX');
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
    
    [v0, g] = f(x);
    step = -g;
    
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

