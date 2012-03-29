function f = tikregf(c)
%TIKREGF Function handle for regularization
%
%   f = TIKREGF(c);
%   f = TIKREGF(Q);
%
%       Constructs and returns a function handle for Tikhonov
%       regularization.
%
%       Tikhonov regularization function is generally defined to be
%
%       f(x) = (1/2) * x' * Q * x.
%
%       Here, Q is a positive semi-definite matrix. Often times, a
%       simplified version, weighted L2-regularization, is used
%
%       f(x) = (1/2) * sum_i c_i (x_i)^2.
%
%       
%       Input arguments:
%       - c:        the coefficients for (weighted) L2 regularization,
%                   which can be in either of the following forms:
%                   - a scalar: all coeffficients are the same
%                   - a vector that explicitly give all coefficients
%
%       - Q:        the coefficient matrix for Tikhonov regularization
%
%       The output f is a function handle, which can be used as follows
%
%           v = f(x);
%           [v, g] = f(x);
%           [v, g, H] = f(x);
%
%       Here, x, the input to f, should be a matrix of size d x n, where
%       n is the number of parameters packed in x (e.g. in multi-class
%       classifier training, n can be greater than 1).
%
%       f returns:
%       - v:        the function value evaluated at x
%       - g:        the gradient vector at x (if needed)
%       - H:        the Hessian matrix at x (if needed)
%

% Created by Dahua Lin, on Jan 14, 2012
%

%% verify input arguments

cty = 0;
if isfloat(c) && ndims(c) == 2
    if isscalar(c)
        cty = 1;        
    elseif isvector(c)
        cty = 2;
        if size(c, 2) > 1
            c = c.';
        end        
    elseif size(c,1) == size(c,2)
        cty = 3; 
        Q = c;
    end    
end

if ~cty
    error('tikregf:invalidarg', 'The input to tikregf is invalid.');
end

%% main

if cty == 1
    f = @(x) reg_L2_sca(x, c);
elseif cty == 2
    f = @(x) reg_L2_vec(x, c);
else
    f = @(x) reg_Q(x, Q);
end

%% regularization functions

function [v, g, H] = reg_L2_sca(x, c)

v = (c/2) * (x' * x);
if nargout >= 2
    g = c * x;
end
if nargout >= 3
    H = diag(c * ones(1, numel(x)));
end



function [v, g, H] = reg_L2_vec(x, c)

g = c .* x;
v = (x' * g) / 2;
if nargout >= 3
    H = diag(c);
end


function [v, g, H] = reg_Q(x, Q)

g = Q * x;
v = (x' * g) / 2;
H = Q;

