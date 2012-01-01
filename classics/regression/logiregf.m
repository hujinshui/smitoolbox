function f = logiregf(X, y, w, rc)
%LOGIREGF Logistic regression objective function
%
%   The objective function of logistic regression is given by
%
%       f(theta) = - \sum_i w_i ( y_i log(p_i) + (1 - y_i) log(1 - p_i) )
%                + (rc / 2) * ||theta||^2
%
%   Here, p_i = 1 / (1 + exp(-z_i)) and z_i = theta' * x_i + theta0
%
%   f = LOGIREGF(X, y);
%   f = LOGIREGF(X, y, w);
%   f = LOGIREGF(X, y, [], rc);
%   f = LOGIREGF(X, y, w, rc);
%       
%       The function returns a function handle f that represents the
%       objective function as formalized above.
%
%       Input arguments:
%       - X:        The sample matrix of size d x n.
%
%       - y:        The indicator vector of length n.
%                   y(i) = 1 indicates that the i-th sample in X is a 
%                   postive sample, and y(i) = 0 indicates that it is
%                   a negative sample.
%
%                   Actually, y(i) can be any real value in [0, 1], in 
%                   which case, it may represent a soft tendency rather
%                   than a hard assignment.                   
%
%       - w:        The weights of the samples. If all samples have the
%                   same weight, then w can be empty or omitted. 
%                   Otherwise, w should be a vector of length n. 
%
%       - rc:       The regularization coefficient.
%
%       The output argument f is a function handle, which can be invoked 
%       as follows:
%
%           v = f(a);       
%           [v, g] = f(a);
%
%       Here, f takes as input the parameter a, in form of [theta; theta0]
%       and returns the objective value v, or optionally the gradient g.
%       This function handle can be used as an objective function in
%       numeric optimization.
%

%   History
%   -------
%       - Created by Dahua Lin, on Jan 25, 2011.
%       - Modified by Dahua Lin, on Jan 1, 2012.
%

%% verify input arguments

if ~(isfloat(X) && ndims(X) == 2)
    error('logiregf:invalidarg', 'X should be a numeric matrix.');
end
n = size(X, 2);

if ~((islogical(y) || isnumeric(y)) && isvector(y) && numel(y) == n)
    error('logiregf:invalidarg', ...
        'y should be a logical or numeric vector of length n');
end
if size(y, 1) > 1
    y = y.';
end
if ~isfloat(y)
    y = double(y);
end

if nargin < 3 || isempty(w)
    w = [];
else
    if ~(isfloat(w) && isvector(w) && numel(w) == n)
        error('logiregf:invalidarg', 'w should be a vector of length n.');
    end
end

if nargin < 4
    rc = 0;
else
    if ~(isfloat(rc) && isscalar(rc) && rc >= 0)
        error('logiregf:invalidarg', 'rc should be a non-negative scalar.');
    end
end


%% main

Xa = [X; ones(1, n)];
f = reglossmin_objfun(Xa, y, 1, w, @logireg_loss, rc, @parreg);


function [v, g] = parreg(a)

theta = a(1:end-1);
v = norm(theta)^2 / 2;
if nargout >= 2
    g = a;
    g(end) = 0;
end

