function f = rbstlr(X, Y, w, rho, s)
% Robust linear regression (M-estimation)
%
%   This function returns an objective function of the following 
%   optimization problem:
%
%       minimize   sum_i rho(||x_i' * a - y_i|| / s) * (s^2)
%
%   Here, rho is a robust loss function that converts the residue to
%   a cost value, s is a scale parameter.
%
%   f = rbstlr(X, Y);
%   f = rbstlr(X, Y, w);
%   f = rbstlr(X, Y, w, rho);
%   f = rbstlr(X, Y, w, rho, s);
%       
%       The function returns a function handle f that represents the
%       optimization objective as formalized above.
%
%       Input arguments:
%       - X:        the design matrix. Suppose there are n input samples,
%                   each with q components, then X should be an matrix
%                   of size n x d.
%
%       - Y:        the response matrix. The size of y should be n x q.
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
%       The output argument f is a function handle, which can be invoked 
%       as follows:
%
%           v = f(a);       
%           [v, g] = f(a);
%           [v, g, H] = f(a);
%
%       Here, f takes as input the parameter a, and returns the objective
%       value v, or optionally the gradient g, and Hessian matrix H.
%       This function handle can be used as an objective function in
%       numeric optimization.
%

%   History
%   -------
%       - Created by Dahua Lin, on Jan 5, 2011
%       - Modified by Dahua Lin, on Jan 6, 2011
%           - supports the case with q > 1.
%       - Modified by Dahua Lin, on Jan 22, 2011
%           - now returns an objective function
%       - Modified by Dahua Lin, on Jan 24, 2011
%           - reimplemented based on genlr
%

%% verify input arguments

if ~(isfloat(X) && ndims(X) == 2)
    error('rbstlr:invalidarg', 'X should be a numeric matrix.');
end
n = size(X, 1);

if ~(isfloat(Y) && ndims(Y) == 2 && size(Y, 1) == n)
    error('rbstlr:invalidarg', 'y should be a numeric matrix with n rows');
end
q = size(Y, 2);

if nargin < 3 || isempty(w)
    w = [];
else
    if ~(isfloat(w) && isvector(w) && numel(w) == n)
        error('rbstlr:invalidarg', 'w should be a vector of length n.');
    end
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
   

%% main

% determine s

if isempty(s)
    s = decide_s(X, Y, w);
end

% form objective function

h = @(u, y) rblr_h(u, y, rho, s);
f = genlr(h, q, X, Y, w);


%% The h-function
    
function [v, D1, D2] = rblr_h(U, Y, rho, s)


[n, q] = size(U);

E = U - Y;
if q == 1
    r2 = E.^2;
    r = abs(E);    
else
    r2 = dot(E, E, 2);
    r = sqrt(r2);
end

rs = r * (1/s);

if nargout <= 1
    v = rho(rs);    
elseif nargout <= 2
    [v, dv1] = rho(rs);
else
    [v, dv1, dv2] = rho(rs);
end

v = v * (s^2);

if nargout >= 2
    ga1 = dv1 ./ rs;
    if q == 1
        D1 = ga1 .* E;
    else
        D1 = bsxfun(@times, ga1, E);
    end
end

if nargout >= 3
    if q == 1
        D2 = dv2; 
    else
        D2 = zeros(n, q, q);
        c2 = (dv2 - ga1) ./ r2;
        
        
        for i = 1 : q
            ei = E(:,i);            
            for j = 1 : i                
                ej = E(:,j);                
                
                h2 = c2 .* (ei .* ej);
                if i == j
                    h2 = h2 + ga1;
                end
                
                D2(:,i,j) = h2;
                if i ~= j
                    D2(:,j,i) = h2;
                end                
            end
        end
    end
end
    
    
%% Auxiliary function

function s = decide_s(X, y, w)

q = size(y, 2);
a0 = llsq(X, y, w);
if q == 1
    mean_dev = mean(abs(X * a0 - y));
else
    e0 = X * a0 - y;
    mean_dev = mean(sqrt(dot(e0, e0, 2)));
end
s = 2 * mean_dev;




    
