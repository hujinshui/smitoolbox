function f = logireg(X, y, w)
% Construct the objective for Logistic regression
%
%   The objective function of logistic regression is given by
%
%       f(theta) = - \sum_i w_i ( y_i log(p_i) + (1 - y_i) log(1 - p_i) )
%
%   Here, p_i = 1 / (1 + exp(-theta' * x_i))
%
%   f = logireg(X, y);
%   f = logireg(X, y, w);
%       
%       The function returns a function handle f that represents the
%       objective function as formalized above.
%
%       Input arguments:
%       - X:        the design matrix. Suppose there are n input samples,
%                   each with q components, then X should be an matrix
%                   of size n x d.
%
%       - y:        the indicator vector of size n x 1. 
%                   y(i) = 1 indicates that the i-th sample in X is a 
%                   postive sample, and y(i) = 0 indicates that it is
%                   a negative sample.
%                   Actually, y(i) can be any real value in [0, 1], in 
%                   which case, it may represent a soft tendency rather
%                   than a hard assignment.                   
%
%       - w:        The weights of the samples. If all samples have the
%                   same weight, then w can be empty or omitted. 
%                   Otherwise, w should be a vector of length n. 
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
%       - Created by Dahua Lin, on Jan 25, 2001.
%

%% verify input arguments

if ~(isfloat(X) && ndims(X) == 2)
    error('logireg:invalidarg', 'X should be a numeric matrix.');
end
n = size(X, 1);

if ~((islogical(y) || isnumeric(y)) && isequal(size(y), [n 1]))
    error('logireg:invalidarg', ...
        'y should be a logical or numeric vector of size n x 1');
end
if ~isfloat(y)
    y = double(y);
end

if nargin < 3 || isempty(w)
    w = [];
else
    if ~(isfloat(w) && isvector(w) && numel(w) == n)
        error('logireg:invalidarg', 'w should be a vector of length n.');
    end
end

%% main

f = genlr(@logir_h, 1, X, y, w);


%% The h-function

function [v, D1, D2] = logir_h(u, y)

pp = 1 ./ (1 + exp(-u));
pn = 1 - pp;

v = - (y .* log(pp) + (1-y) .* log(pn));

if nargout >= 2
    D1 = pp - y;
end

if nargout >= 3
    D2 = pp .* pn;
end

