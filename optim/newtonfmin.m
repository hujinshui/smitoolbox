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
%       - Display:  Level of display. {'off'}|'iter'|'final'|'notify'
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
%

%% verify input 

verify_optim_options('newtonfmin', options, 'MaxIter', 'TolFun', 'TolX', 'Display');
disp_level = get_display_level('newtonfmin', options);

%% main

x = x0;
converged = false;
it = 0;

while ~converged && it < options.MaxIter
    
    it = it + 1;
   
    [v0, g, H] = f(x);
    step = - H \ g;
    
    [x, v, dx] = linesearch(f, x0, v0, step, beta, minstep);
    
    ch = v - v0;
    nrm_dx = norm(dx);
    converged = abs(ch) < options.TolFun && nrm_dx < options.TolX;  
    
    if disp_level >= 3
        fprintf('Iter %d: objv: %.4g => %.4g [ch = %.4g, norm(dx) = %.4g]', ...
            it, v0, v, ch, nrm_dx);        
    end
end

if disp_level >= 1
    if converged && disp_level >= 2
        fprintf('The optimization converges with %d iterations.\n', it);
    elseif ~converged
        fprintf('The optimization terminates without convergence at iteration %d.\n', it);
    end
end


if nargout >= 2
    info = struct( ...
        'FunValue', v, ...
        'LastChange', ch, ...
        'LastMove', dx, ...
        'IsConverged', converged, ...
        'NumIters', it);
end



