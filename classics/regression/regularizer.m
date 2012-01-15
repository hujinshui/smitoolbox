function f = regularizer(type, coeffs)
%REGULARIZER Function handle for regularization
%
%   f = REGULARIZER(type, coeffs);
%
%       Constructs and returns a regularization function handle of the
%       specified type.
%       
%       Input arguments:
%       - type:     The type of regularizer, which can be either of 
%                   the following strings
%                   - 'L1':     The L1-norm regularizer
%                               f(x) = sum_i c_i * |x_i|
%
%                   - 'L2':     The L2-norm regularizer
%                               f(x) = (1/2) * sum_i c_i * (x_i)^2
%
%                   - 'Q':      The quadratic regularizer
%                               f(x0 = (1/2) * x' * Q * x.
%
%       - coeffs:   The regularization coefficients. 
%                   If type is 'L1' or 'L2', it can be in either
%                   of the following forms:
%
%                   - a scalar: all coefficients have the same values
%                               as of the scalar.
%
%                   - a vector of length d. Here, d is the dimension of
%                     the parameters space.
%
%                   - a cell array in form of {c, c0}. Here, in the 
%                     solution x, x(1:end-1,:) represents the coefficient
%                     vector(s), and x(end,:) represents the offsets.
%                     And, c is the coefficient for the former part, while
%                     c0 is for the latter.
%
%                   If type is 'Q', it can be in either of the following
%                   forms:
%                   - Q, the d x d coefficient matrix. 
%                   - A cell array in form of {Q, c0}, where Q is used
%                     to regularize x(1:end-1,:).
%
%       The output f is a function handle, which can be used as follows
%
%           v = f(x);
%           [v, g] = f(x);
%
%       Here, x, the input to f, should be a matrix of size d x n, where
%       n is the number of parameters packed in x (e.g. in multi-class
%       classifier training, n can be greater than 1).
%
%       f returns v the regularization term value for x (if n > 1, then
%       v is the sum of the regularization values for all parameters).
%       f also returns g, the gradient at x, when there are two output
%       arguments.
%

% Created by Dahua Lin, on Jan 14, 2012
%

%% verify input arguments

if strcmpi(type, 'L1')
    ty = 1;
elseif strcmpi(type, 'L2')
    ty = 2;
elseif strcmpi(type, 'Q')
    ty = 3;
else
    error('regularizer:invalidarg', 'The regularizer type is invalid.');
end

if ty == 1 || ty == 2
    if isnumeric(coeffs)
        if ~(isfloat(coeffs) && isvector(coeffs))
            error('regularizer:invalidarg', ...
                'The coeffs should be a real scalar or vector.');
        end
        if isscalar(coeffs)
            c = coeffs;
        else
            c = coeffs;
            if size(c, 2) > 1
                c = c.';
            end
        end
        
    elseif iscell(coeffs) && numel(coeffs) == 2
        c = coeffs{1};
        c0 = coeffs{2};
        
        if ~(isfloat(c) && isscalar(c) && isfloat(c0) && isscalar(c0))
            error('regularizer:invalidarg', 'The regularizer coeffs is invalid.');
        end
        
        if size(c, 2) > 1
            c = c.';
        end
        c = [c; c0];
    else
        error('regularizer:invalidarg', 'The regularizer coeffs is invalid.');
    end
else
    if isnumeric(coeffs)
        Q = coeffs;        
    elseif iscell(coeffs) && numel(coeffs) == 2
        Q = coeffs{1};
        c0 = coeffs{2};
        
        Q = blkdiag(Q, c0);
    end
    
    if ~(isfloat(Q) && ndims(Q) == 2 && size(Q,1) == size(Q,2))
        error('regularizer:invalidarg', ...
            'Q should be a real square matrix.');
    end
end
        
%% main

if ty == 1
    f = @(x) reg_L1(x, c);
elseif ty == 2
    f = @(x) reg_L2(x, c);
else
    f = @(x) reg_Q(x, Q);
end


%% regularization functions

% L1

function [v, g] = reg_L1(x, c)

v = calc_v(c, abs(x));
if nargout >= 2
    g = calc_g(c, sign(x));
end


% L2

function [v, g] = reg_L2(x, c)

v = calc_v((c/2), x.^2);
if nargout >= 2
    g = calc_g(c, x);
end

% Q

function [v, g] = reg_Q(x, Q)

g = Q * x;
v = 0.5 * sum(sum(x .* g, 1));



%% auxiliary functions

function v = calc_v(c, u)

if size(u, 2) == 1
    if isscalar(c)
        v = c * sum(u);
    else
        v = sum(c .* u);
    end
else
    if isscalar(c)
        v = c * sum(sum(u, 1));
    else
        v = sum(sum(bsxfun(@times, c, u), 1));
    end
end


function g = calc_g(c, u)

if isscalar(c) || size(u, 2) == 1
    g = c .* u;
else
    g = bsxfun(@times, c, u);
end
    
