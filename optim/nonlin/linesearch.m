function [x, v, dx, fcnt] = linesearch(f, x0, v0, step, beta, minstep)
% Performs line search as a step in numerical optimization
%
%   [x, v, dx, fcnt] = linesearch(f, x, step, beta, minstep);
%       performs line search along the direction of step to locate the
%       minima.
%
%       Inputs:
%       - f:        the objective function to be minimized
%       - x0:       the location from which the search starts
%       - v0:       the objective value at x0: f(x0)
%       - step:     the initial moving step
%       - beta:     the exponential decay coefficient
%       - minstep:  the minimum allowable L2-norm of a step
%
%       Outputs:
%       - x:        the located point (x0 + dx)
%       - v:        the objective value at located point (x)
%       - dx:       the actual step 
%       - fcnt:     the number of times that f is invoked
%
%       The search procedure runs as follows. It first take dx = step,
%       if f(x + dx) < f(x) then it is done, otherwise, it decreases dx
%       as dx = dx * beta. This process continues until the condition
%       f(x + dx) < f(x) is met, or ||dx|| < minstep. In latter case,
%       we simply return with dx set to a zero vector.
%
%   Remarks
%   --------
%       - This function is designed to be used by an optimization
%         algorithm.
%       - As this function is on a critical path of a typical optimization
%         method, no argument checking is performed for the sake of
%         efficiency. 
%           
%       

% Created by Dahua Lin, on Aug 2, 2010
%

%% main

dx = step;
x = x0 + dx;
v = f(x);
fcnt = 1;

while v >= v0 && norm(dx) > minstep    
    dx = dx * beta;    
    x = x0 + dx;
    v = f(x);
    fcnt = fcnt + 1;
end

if v >= v0
    dx = 0;
    x = x0;
    v = v0;
end


