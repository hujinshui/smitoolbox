function objf = reglossmin_objfun(X, Y, K, w, lossfun, rc, regfun)
%REGLOSSMIN_OBJFUN Regularized Loss Minimization Objective function
%
%   Generally, a regularized loss minimization problem is formulated as
%
%       minimize sum_{i=1}^n w_i loss(theta' * x_i, y_i) + rc * reg(theta).
%
%   Suppose the sample space dimension is d, and the response space
%   dimension is q. Then each x_i is a d x 1 vector, and each y_i is 
%   a q x 1 vector. The parameter theta can be a d x 1 vector or a 
%   d x K matrix.
%
%
%   objf = reglossmin_objfun(X, Y, K, w, lossfun);
%   objf = reglossmin_objfun(X, Y, K, w, lossfun, rc);
%   objf = reglossmin_objfun(X, Y, K, w, lossfun, rc, regfun);
%
%       Constructs the objective function for regularized loss
%       minimization.
%
%       Suppose there are n training pair of samples and responses.
%
%       Input arguments:
%       - X:        the sample (feature) matrix, size: d x n.
%
%       - Y:        the response matrix, size: q x n.
%
%       - K:        the number of columns in the parameter matrix.
%
%       - w:        the sample weights, size: 1 x n.
%                   w can also be [], indicating all sample weights are 1.
%
%       - lossfun:  the function handle of the loss function. It will be
%                   invoked in the following way:
%
%                   [v, g] = lossfun(theta' * X, Y).
%
%                   Here, v should be a 1 x n row vector containing the 
%                   objective values, and g should be a K x n matrix
%                   containing the gradient w.r.t. the first argument.
%
%       - rc:       the regularization term coefficient.
%
%       - regfun:   the regularization function, which will be invoked 
%                   in the following way:
%
%                   [v, g] = regfun(theta).
%
%                   Here, v should be a scalar representing the objective
%                   value, and g should be a d x K matrix, the same size 
%                   as the parameter theta.
%
%                   If regfun is omitted, L2 regularization is used, as
%
%                       regfun(theta) = (1/2) * ||theta||^2.
%
%                   
%       Output arguments:
%       - objf:     an objective function handle which can be used as
%                   follows:
%
%                       v = objf(theta).
%                       [v, g] = objf(theta).
%
%                   Here, v is the objective value, and g is the gradient
%                   vector of length d x K.
%
%       Note that objf can serve as an input to optimization functions
%       such as fminunc. The solution will be a vector of length d x K,
%       which is a concatenation of all columns of theta.
%       

%   History
%   -------
%       - Created by Dahua Lin, on Jan 24, 2011
%

%% verify inputs

if ~(isfloat(X) && isreal(X) && ndims(X) == 2)
    error('reglossmin_objfun:invalidarg', 'X should be a real matrix.');
end

if ~(isfloat(Y) && isreal(Y) && ndims(Y) == 2)
    error('reglossmin_objfun:invalidarg', 'Y should be a real matrix.');
end

if ~(isnumeric(K) && isscalar(K) && K == fix(K) && K >= 1)
    error('reglossmin_objfun:invalidarg', 'K should be a positive integer.');
end

if ~isempty(w)
    if ~(isfloat(w) && isreal(w) && isvector(w) && length(w) == n)
        error('reglossmin_objfun:invalidarg', ...
            'w should be a real vector of length n.');
    end
    if size(w, 2) > 1; w = w.'; end  % turn into column
end

if ~isa(lossfun, 'function_handle')
    error('reglossmin_objfun:invalidarg', ...
        'lossfun should be a function handle.');
end

if nargin >= 6
    if ~(isfloat(rc) && isreal(rc) && isscalar(rc) && rc >= 0)
        error('reglossmin_objfun:invalidarg', ...
            'rc should be a non-negative real scalar.');
    end    
    
    if nargin >= 7
        if ~isa(regfun, 'function_handle')
            error('reglossmin_objfun:invalidarg', ...
                'regfun should be a function handle.');            
        end
    else
        regfun = @regL2;
    end
else
    rc = 0;
    regfun = [];
end

%% main

objf = @(theta) rlmin_obj(theta, X, Y, K, w, lossfun, rc, regfun);

%% The objective function

function [v, g] = rlmin_obj(theta, X, Y, K, w, lossfun, rc, regfun)

d = size(X, 1);
if K > 1
    theta = reshape(theta, d, K);
end

Z = theta' * X;

out_g = nargout >= 2;

if out_g
    [L, Gz] = lossfun(Z, Y);  % Gz: K x n
else
    L = lossfun(Z, Y);
end

if isempty(w)
    t_lv = sum(L);
    if out_g
        t_lg = X * Gz';
    end
else
    t_lv = L * w;
    if out_g
        if K == 1
            t_lg = X * (Gz' .* w);
        else
            t_lg = X * bsxfun(@times, Gz', w);
        end
    end    
end

if out_g && K > 1
    t_lg = t_lg(:);
end

if rc > 0
    if out_g
        [rv, rg] = regfun(theta);    
        if size(rg, 2) ~= 1
            rg = rg(:);
        end
    else
        rv = regfun(theta);
    end
    
    v = t_lv + rc * rv;    
    if out_g
        g = t_lg + rc * rg;
    end
else
    v = t_lv;    
    if out_g
        g = t_lg;
    end
end



%% The regularization function

function [v, g] = regL2(a)

if size(a, 2) > 1
    a = a(:);
end

v = norm(a)^2 / 2;
g = a;


