function f = comb_objfun(varargin)
%COMB_OBJFUN Linear combination of objective functions
%
%   An objective function f is a function that support the following use:
%
%   (1) v = f(x);
%   (2) [v, g] = f(x);  (optional)
%   (3) [v, g, H] = f(x);  (optional)
%
%   Here, x is a feasible solution, and v, g, H are respective the 
%   function value, gradient, and Hessian matrix evaluated at x.
%
%
%   f = comb_objfun(c1, f1);
%       returns an objective function f, which is equivalent to c1 * f1.
%
%   f = comb_objfun(c1, f1, c2, f2);
%       returns an objective function f, as c1 * f1 + c2 * f2
%
%   f = comb_objfun(c1, f1, c2, f2, c3, f3, ...);
%       returns an objective function f, as a linear combination of
%       multiple objective functions.
%

%% verify input arguments

if isempty(varargin)
    error('comb_objfun:invalidarg', 'There is no input to comb_objfun');
end

coeffs = varargin(1:2:end);
funs = varargin(2:2:end);

n = numel(coeffs);
if ~(numel(funs) == n && ...
        all(cellfun(@(x) isfloat(x) && isscalar(x), coeffs)) && ...
        all(cellfun(@(f) isa(f, 'function_handle'), funs)) )
    error('comb_objfun:invalidarg', 'The inputs to comb_objfun are invalid.');
end


%% main

if n == 1
    c1 = coeffs{1};
    f1 = funs{1};
    f = @(x) comb_objfun1(x, c1, f1);
    
elseif n == 2
    c1 = coeffs{1};
    c2 = coeffs{2};
    f1 = funs{1};
    f2 = funs{2};
    f = @(x) comb_objfun2(x, c1, f1, c2, f2);
    
else 
    c = vertcat(coeffs{:});
    f = @(x) comb_objfun_multi(x, c, funs);
    
end


%% objective functions

function [v, g, H] = comb_objfun1(x, c1, f1)

nout = nargout;
if nout <= 1
    v = c1 * f1(x);
elseif nout == 2
    [v1, g1] = f1(x);
    v = c1 * v1;
    g = c1 * g1;
elseif nout == 3
    [v1, g1, H1] = f1(x);
    v = c1 * v1;
    g = c1 * g1;
    H = c1 * H1;
end

function [v, g, H] = comb_objfun2(x, c1, f1, c2, f2)

nout = nargout;
if nout <= 1
    v = c1 * f1(x) + c2 * f2(x);
elseif nout == 2
    [v1, g1] = f1(x);
    [v2, g2] = f2(x);
    v = c1 * v1 + c2 * v2;
    g = c1 * g1 + c2 * g2;
elseif nout == 3
    [v1, g1, H1] = f1(x);
    [v2, g2, H2] = f2(x);
    v = c1 * v1 + c2 * v2;
    g = c1 * g1 + c2 * g2;
    H = c1 * H1 + c2 * H2;
end

function [v, g, H] = comb_objfun_multi(x, c, funs)

nout = nargout; 
m = numel(funs);

if nout <= 1
    f1 = funs{1};
    v = c(1) * f1(x);
    for k = 2 : m
        fk = funs{k};
        v = v + c(k) * fk(x);
    end
    
elseif nout == 2
    f1 = funs{1};
    [v1, g1] = f1(x);
    v = c(1) * v1;
    g = c(1) * g1;
    for k = 2 : m
        fk = funs{k};
        [vk, gk] = fk(x);
        v = v + c(k) * vk;
        g = g + c(k) * gk;
    end
    
elseif nout == 3
    f1 = funs{1};
    [v1, g1, H1] = f1(x);
    v = c(1) * v1;
    g = c(1) * g1;
    H = c(1) * H1;
    for k = 2 : m
        fk = funs{k};
        [vk, gk, Hk] = fk(x);
        v = v + c(k) * vk;
        g = g + c(k) * gk;
        H = H + c(k) * Hk;
    end
    
end




