function f = genlr(h, q, X, Y, w)
% Construct a generalized linear regression objective function
%
%   A generalized linear regression problem is to minimize the following
%   objective function:
%
%       f(theta) = sum_i w_i h(theta'*x_i, y_i).
%
%   Here, h is a function that measures the compatibility between the
%   linear predictor and the response y.
%
%   f = genlr(h, q, X, Y);
%   f = genlr(h, q, X, Y, w);
%
%       This statement returns the function handle f that captures the
%       objective function as formalized above. 
%
%       Input arguments:
%
%       - h:    The compatibility function for each term. It is a function
%               handle that can be called as follows
%
%                   v = h(U, Y);
%
%               Here, h takes as input the linear predictors in U, a matrix
%               of size n x q, and the response matrix Y, and returns the 
%               function values in an n x 1 vector v.
%              
%               If the objective is to be used in numeric optimization 
%               that would use gradient or Hessian. Then it should
%               support the following way of function call:
%
%                   [v, D1] = h(U, Y);
%                   [v, D1, D2] = h(U, Y);
%
%               Here, D1 is a matrix of size n x q, each row of D1
%               gives the first-order derivative with respect to u for
%               a sample. D2 is an array of size n x 1 or n x q x q,
%               where D2(i,:,:) gives the second-order derivative with
%               respect to u for a sample. 
%
%               Note that h could be implemented in a vectorized way to 
%               make the computation more efficient.
%
%       - q:    The number of columns in parameters. The size of the
%               parameter theta is d x q.
%
%       - X:    The design matrix. Suppose there are n samples, and each
%               sample has d independent variables (explanatory variables).
%               Then X should be a matrix of size n x d, with X(i,:) 
%               corresponding to the i-th sample.
%
%       - Y:    The response matrix. Suppose there are q dependent 
%               variables for each sample. Then Y should be a matrix with 
%               n rows, and Y(i,:) corresponds to the i-th sample.
%
%       - w:    The weights of terms. If it is omitted or left empty, it
%               indicates that each term has a weight 1. It can also be
%               a vector of length n.
%
%       Output argument:
%
%       - f:    The function handle that represents the objective. It 
%               can be used in the following way:
%
%                   v = f(theta);
%                   [v, g] = f(theta);
%                   [v, g, H] = f(theta);
%
%               Here, v, g, and H are respectively the function value,
%               the gradient, and the Hessian matrix of f evaluated at 
%               theta. Note that to compute g and H, the supplied function
%               h should support the computation of 1st and 2nd order
%               derivatives.
%
%   Remarks
%   -------
%       - One can use an unconstrained nonlinear optimization function
%         to solve the optimal parameter. For example, you can write
%
%           f = genlr( ... );
%           theta = fminunc(f, theta0(:));
%           theta = reshape(a, [d q]);
%
%         Here, fminunc is used together with genlr to solve the optimal
%         parameter theta.
%

%   History
%   -------
%       - Created by Dahua Lin, on Jan 24, 2011
%

%% verify input arguments

if ~isa(h, 'function_handle')
    error('genlr:invalidarg', 'h should be a function handle.');
end

if ~(isnumeric(q) && isscalar(q) && q == fix(q) && q >= 1)
    error('genlr:invalidarg', 'q should be a positive integer scalar.');
end

if ~(isfloat(X) && ndims(X) == 2 && isreal(X))
    error('genlr:invalidarg', 'X should be a real matrix.');
end
n = size(X, 1);

if ~(isnumeric(Y) && ndims(Y) == 2)
    error('genlr:invalidarg', 'Y should be a numeric matrix.');
end

if nargin < 4 || isempty(w)
    w = [];
else
    if ~(isfloat(w) && isvector(w) && numel(w) == n)
        error('genlr:invalidarg', ...
            'w should be either empty or a vector of length n.');
    end
    if size(w, 2) > 1; w = w.'; end     % turn into a column        
end

%% main

f = @(theta) glr(theta, h, q, X, Y, w);


%% The objective function

function [v, g, H] = glr(theta, h, q, X, Y, w)

d = size(X, 2);
if ~isequal(size(theta), [d, q])
    theta = reshape(theta, d, q);
end

U = X * theta;  % linear predictor

% per-term evaluation

if nargout <= 1
    V = h(U, Y);
elseif nargout <= 2
    [V, D1] = h(U, Y);
else
    [V, D1, D2] = h(U, Y);
end

% integration

% function value

if isempty(w)
    v = sum(V);
else
    v = w' * V;
end

% gradient

if nargout >= 2
    if isempty(w)
        g = X' * D1;
    else
        g = X' * bsxfun(@times, w, D1);
    end
    
    if q > 1
        g = g(:);
    end
end

% Hessian

if nargout >= 3
    if q == 1
        H = glr_make_H(X, w, D2);
        
    else % q > 1
        H = zeros(d*q, d*q);
        
        for i = 1 : q
            
            si = (i-1)*d + 1;
            ei = i*d;            
            
            for j = 1 : i
                
                sj = (j-1)*d + 1;
                ej = j*d;
                                
                H_ij = glr_make_H(X, w, D2(:,i,j));
                
                H(si:ei, sj:ej) = H_ij;
                if j < i
                    H(sj:ej, si:ei) = H_ij';
                end
            end
        end        
    end
end


%% Auxiliary function

function H = glr_make_H(X, w, D2)

if isempty(w)
    w2 = D2;
else
    w2 = w .* D2;
end

H = X' * bsxfun(@times, w2, X);
H = 0.5 * (H + H');
    
