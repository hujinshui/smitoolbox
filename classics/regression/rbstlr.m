function [a, info] = rbstlr(X, y, f, s, varargin)
% Robust linear regression (M-estimation)
%
%   This function solves the following problem
%
%       minimize sum_i rho(||x_i' * a - y_i|| / s) * (s^2)
%
%   Here, rho is a robust loss function that converts the residue to
%   a cost value, s is a scale parameter.
%
%   a = rbstlr(X, y, f);
%   a = rbstlr(X, y, f, s, ...);
%   [a, info] = rbstlr( ... );
%   
%       solves te problem as formalized above. Suppose there are n 
%       independent variables (x_1, ..., x_n), each with d components.
%       Then X is an n x d matrix, with X(i, :) corresponding to x_i.
%       And, y should be an n x 1 vector.
%
%       The robust error function is specified by f. Here, f can be
%       a string giving the name of a pre-defined loss function.
%       (see the help of robustloss for a list of available names).
%
%       One can also plug in a user-specified loss function, in form of
%       a function handle. Here, f should support the following usages:
%
%           % returns function value, 1st, and 2nd order derivatives
%           [v, dv1, dv2] = f(z);  
%
%           % returns the weight value, w = dv1 / v
%           w = f(z, 'w');
%
%       It takes a matrix z, and returns the results in matrices of the 
%       the same size. 
%
%       The 4th argument s is the scale parameter. If s is omitted, it 
%       will be determined as 2 * average residue. 
%
%       In addition to these parameters, one can specify further options
%       in name/value pairs to control the optimization process:
%
%       - 'Method':     which method to use. The value can be:
%                       - 'newton': Use newton's method (default)
%                       - 'irls':   Iterative reweighted least square
%       - 'MaxIter':    maximum number of iterations {100}
%       - 'TolFun':     the tolerance of objective function at convergence
%                       {1e-10}. 
%       - 'TolX':       the tolerance of solution change at convergence
%                       {1e-8}
%       - 'Weights':    the sample weights {1}.
%       - 'Init':       initial guess of the coefficient vector a
%                       If not explicitly given, the function initializes
%                       a by ordinary linear regression.
%       - 'L2R':        L2 regularization coefficient {0}
%       - 'Display':    the level of displaying: {'none'}|'proc'|'iter'
%
%       In the output arguments, a is the resultant coefficient vector.
%       The second argument info is a struct that contains the information
%       about the optimization procedure.
%

%   History
%   -------
%       - Created by Dahua Lin, on Jan 5, 2011
%

%% verify input arguments

error(nargchk(3, inf, nargin));

if ~(isfloat(X) && ndims(X) == 2)
    error('rbstlr:invalidarg', 'X should be a numeric matrix.');
end
[n, d] = size(X);

if ~(isfloat(y) && isequal(size(y), [n 1]))
    error('rbstlr:invalidarg', 'y should be a numeric vector of size n x 1.');
end

if ischar(f)
    f = robustloss(f);
else
    if ~isa(f, 'function_handle')
        error('rbstlr:invalidarg', ...
            'f should be either a string or a function handle.');
    end
end

if nargin < 4
    s = [];
else
    if ~(isempty(s) || (isfloat(s) && isscalar(s) && s > 0))
        error('rbstlr:invalidarg', 's should be either empty or a positive scalar.');
    end
end

[mcode, opts, w, a0, r2] = parse_options(n, d, varargin);
   

%% main

% Initialize

if isempty(a0)
    a0 = llsq(X, y, w, r2);
end

if isempty(s)
    s = 2 * mean(abs(X * a0 - y));
end

% delegate to solve

if mcode == 1   % use newton
    if nargout < 2
        a = rlr_newton(X, y, w, r2, f, s, opts, a0);
    else
        [a, info] = rlr_newton(X, y, w, r2, f, s, opts, a0);
    end
else            % use irls
    if nargout < 2
        a = rlr_irls(X, y, w, r2, f, s, opts, a0);
    else
        [a, info] = rlr_irls(X, y, w, r2, f, s, opts, a0);
    end
end


%% Core function

function [v, g, H] = objfun(a, X, y, w, r2, f, s)
% objective function

e = (X * a - y) * (1/s);

[V, D1, D2] = f(e);

if isequal(w, 1)
    v = (s^2) * sum(V);
    g = s * (X' * D1);
    H = X' * bsxfun(@times, D2, X);
else
    v = (s^2) * (V' * w);
    g = s * (X' * (D1 .* w));
    H = X' * bsxfun(@times, D2 .* w, X);
end
H = 0.5 * (H + H');

if r2 > 0
    v = v + (0.5 * r2) * (a' * a);
    g = g + r2 * a;
    H = adddiag(H, r2);
end
    

function [a, info] = rlr_newton(X, y, w, r2, f, s, opts, a0)
% solve using newton's method

F = @(c) objfun(c, X, y, w, r2, f, s);

if nargout < 2
    a = newtonfmin(F, a0, opts{:});
else
    [a, info] = newtonfmin(F, a0, opts{:});
end


function [a, info] = rlr_irls(X, y, w, r2, f, s, opts, a0)
% solve using iterative reweighted least square

if isequal(w, 1)
    wf = @(e) f(e, 'w');
else
    wf = @(e) w .* f(e, 'w');
end

if nargout < 2
    a = irls(X, y, wf, s, r2, a0, opts{:});
else
    [a, rw, info] = irls(X, y, wf, s, r2, a0, opts{:}); %#ok<ASGLU>
end



%% option parsing

function [mcode, opts, w, a0, r2] = parse_options(n, d, nvs)

mcode = 1;  % 1 - newton, 2 - irls
w = 1;
a0 = [];
r2 = 0;

maxiter = 100;
tolfun = 1e-10;
tolx = 1e-8;
display = 0;

if ~isempty(nvs)
    ns = nvs(1:2:end);
    vs = nvs(2:2:end);
    if ~(numel(ns) == numel(vs) && iscellstr(ns))
        error('rbstlr:invalidarg', 'The option list is invalid.');
    end
    for i = 1 : numel(ns)
        nam = ns{i};
        v = vs{i};
        switch lower(nam)
            case 'method'
                if ~(ischar(v))
                    error('rbstlr:invalidarg', 'Method should be a string.');
                end
                if strcmpi(v, 'newton')
                    mcode = 1;
                elseif strcmpi(v, 'irls')
                    mcode = 2;
                else
                    error('rbstlr:invalidarg', 'The Method is invalid.');
                end
                
            case 'maxiter'
                if ~(isnumeric(v) && isscalar(v) && v >= 1)
                    error('rbstlr:invalidarg', ...
                        'MaxIter should be a positive integer scalar.');
                end
                maxiter = v;
                
            case 'tolfun'
                if ~(isfloat(v) && isscalar(v) && v > 0)
                    error('rbstlr:invalidarg', ...
                        'TolFun should be a positive scalar.');
                end
                tolfun = v;
                
            case 'tolx'
                if ~(isfloat(v) && isscalar(v) && v > 0)
                    error('rbstlr:invalidarg', ...
                        'TolX should be a positive scalar.');
                end
                tolx = v;
                
            case 'weights'
                if ~(isfloat(v) && isvector(v) && numel(v) == n)
                    error('rbstlr:invalidarg', ...
                        'Weights should be a vector of length n.');
                end
                w = v;
                if size(w, 2) > 1; w = w.'; end              
                
            case 'init'
                if ~(isfloat(v) && isequal(size(v), [d 1]))
                    error('rbstlr:invalidarg', ...
                        'Init should be a vector of size d x 1.');
                end
                a0 = v;
                
            case 'l2r'
                if ~(isfloat(v) && isscalar(v) && v >= 0)
                    error('rbstlr:invalidarg', ...
                        'L2R should be a positive scalar.');
                end
                r2 = v;
                
            case 'display'
                display = v;
                
            otherwise
                error('rbstlr:invalidarg', ...
                    'Invalid option name %s', nam);
        end
    end    
end

if display == 0
    opts = {'MaxIter', maxiter, 'TolFun', tolfun, 'TolX', tolx};
else
    opts = {'MaxIter', maxiter, 'TolFun', tolfun, 'TolX', tolx, 'Display', display};
end




    
