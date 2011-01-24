function f = rbstlr(X, y, w, rho, s, r2)
% Robust linear regression (M-estimation)
%
%   This function returns an objective function of the following 
%   optimization problem:
%
%       minimize   sum_i rho(||x_i' * a - y_i|| / s) * (s^2)
%                + (r2/2) * ||a||^2
%
%   Here, rho is a robust loss function that converts the residue to
%   a cost value, s is a scale parameter.
%
%   f = rbstlr(X, y);
%   f = rbstlr(X, y, w);
%   f = rbstlr(X, y, w, rho);
%   f = rbstlr(X, y, w, rho, s);
%   f = rbstlr(X, y, w, rho, s, r2);
%       
%       The function returns a function handle f that represents the
%       optimization objective as formalized above.
%
%       Input arguments:
%       - X:        the design matrix. Suppose there are n input samples,
%                   each with q components, then X should be an matrix
%                   of size n x d.
%
%       - y:        the response matrix. The size of y should be n x q.
%                   Here, q is the dimension of the response space.
%
%       - w:        The weights of the samples. If all samples have the
%                   same weight, then w can be empty or omitted. 
%                   Otherwise, w should be a vector of length n.       
%
%       - rho:      A robust cost function that converts deviation norm
%                   to a cost value. It can be a name of pre-defined
%                   function (refer to the help of robustloss), or 
%                   a function handle.
%                   If omitted, it is set to the default value 'bisquare'.
%
%       - s:        the scale parameter s. If omitted, it is determined
%                   as twice the mean deviation based on the coefficient
%                   estimated by linear least square.
%
%       - r2:       the L2 regularization coefficient. If omitted, it
%                   is set to zero.
%
%       The output argument f is a function handle, which can be invoked 
%       as follows:
%
%           v = f(a);       
%           [v, g] = f(a);
%           [v, g, H] = f(a);
%
%       Here, f takes as input the parameter a, and returns the objective
%       value v, or optionally the gradient g, and Hessian matrix H.
%       
%       One can use an unconstrained nonlinear optimization function
%       to solve the optimal parameter. For example, you can write
%
%           a = fminunc(rbstlr(X, y), a0(:));
%           a = reshape(a, [d q]);
%
%       it solve the optimal parameter a, give the data X and y and initial
%       guess a0.
%
%   Remarks
%   -------
%       - Note that the output function handle f can produce the Hessian
%         matrix H only when q == 1.
%

%   History
%   -------
%       - Created by Dahua Lin, on Jan 5, 2011
%       - Modified by Dahua Lin, on Jan 6, 2011
%           - supports the case with q > 1.
%       - Modified by Dahua Lin, on Jan 22, 2011
%           - now returns an objective function
%

%% verify input arguments

error(nargchk(2, inf, nargin));

if ~(isfloat(X) && ndims(X) == 2)
    error('rbstlr:invalidarg', 'X should be a numeric matrix.');
end
n = size(X, 1);

if ~(isfloat(y) && ndims(y) == 2 && size(y, 1) == n)
    error('rbstlr:invalidarg', 'y should be a numeric matrix with n rows');
end

if nargin < 3 || isempty(w)
    w = 1;
else
    if ~(isfloat(w) && isvector(w) && numel(w) == n)
        error('rbstlr:invalidarg', 'w should be a vector of length n.');
    end
    if size(w, 2) > 1; w = w.'; end  % turn into a column vector
end

if nargin < 4
    rho = 'bisquare';
end
if ischar(rho)
    rho = robustloss(rho);
else
    if ~isa(rho, 'function_handle')
        error('rbstlr:invalidarg', ...
            'rho should be either a string or a function handle.');
    end
end

if nargin < 5
    s = [];
else
    if ~(isempty(s) || (isfloat(s) && isscalar(s) && s > 0))
        error('rbstlr:invalidarg', 's should be either empty or a positive scalar.');
    end
end

if nargin < 6
    r2 = 0;
else
    if ~(isfloat(r2) && isscalar(r2) && r2 >= 0)
        error('rbstlr:invalidarg', 'r should be a non-negative scalar.');
    end
end
   

%% main

% determine s

if isempty(s)
    s = decide_s(X, y, w, r2);
end

% form objective function

f = @(c) objfun(c, X, y, w, r2, rho, s);



%% Core function

function [v, g, H] = objfun(a, X, y, w, r2, f, s)
% objective function

u_g = nargout >= 2;
u_h = nargout >= 3;
d = size(X, 2);
q = size(y, 2);
if q > 1
    a = reshape(a, d, q);
end

e = (X * a - y) * (1/s);

if q == 1
    [V, D1, D2] = f(e);    
else
    enrm = sqrt(dot(e, e, 2));
    [V, D1] = f(enrm);
    rw = D1 ./ enrm;    
end   
    
if isequal(w, 1)
    v = (s^2) * sum(V);
    if u_g
        if q == 1
            g = s * (X' * D1);
        else
            g = s * (X' * bsxfun(@times, rw, e));
        end
    end
    if u_h
        if q > 1
            error('rbstlr:rterror', 'Cannot compute H when q > 1.');
        end
        H = X' * bsxfun(@times, D2, X);
        H = 0.5 * (H + H');
    end
else
    v = (s^2) * (V' * w);
    if u_g
        if q == 1
            g = s * (X' * (D1 .* w));
        else
            g = s * (X' * bsxfun(@times, rw .* w, e));
        end
    end
    if u_h
        if q > 1
            error('rbstlr:rterror', 'Cannot compute H when q > 1.');
        end
        H = X' * bsxfun(@times, D2 .* w, X);
        H = 0.5 * (H + H');
    end
end

if r2 > 0
    v = v + (0.5 * r2) * (a' * a);
    if u_g
        g = g + r2 * a;
    end
    if u_h
        H = adddiag(H, r2);
    end
end

if u_g && q > 1
    g = g(:);
end


%% Auxiliary function

function s = decide_s(X, y, w, r2)

q = size(y, 2);
a0 = llsq(X, y, w, r2);
if q == 1
    mean_dev = mean(abs(X * a0 - y));
else
    e0 = X * a0 - y;
    mean_dev = mean(sqrt(dot(e0, e0, 2)));
end
s = 2 * mean_dev;




    
